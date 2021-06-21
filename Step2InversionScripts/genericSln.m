classdef genericSln < handle
    %{
6/16/21 Generic Solution
For use with the MCMC algorithm script. This object contains all the
information for a proposed subsurface structure ('solution')
containing some number of layers, each with a thickness (rather, a depth to
the top of the layer) and a resistivity (rho). It also contains all
information for updating itself and checking itself within bounds when a
step is taken in the markov chain.
    Set access to private, properties should only be allowed to be changed 
within this object. If info needs to be retrieved externally, create a
get function for it. If something needs to be set externally, create a set 
function. Get and set functions are at the end of the methods section.
    %}
    properties (Access = private)
        %%%%%%%% Solution physical properties %%%%%%%%%%%%%%%%%%%%%%%%
        numLayers;  % number of layers in proposed solution
        lDepths;  % log-depths (log-meters) to each interface
        lRhos;    % log-resistivities (log-ohm-meters) for each layer
        %%%%%%%%%%%%%%%%%% Solution math properties %%%%%%%%%%%%%%%%%%%%
        var;        % variance \sigma^2, NOT LOG
        residual;   % residual vector, ie misfit from data, ohm-meters
        misfit;     % norm of residual ohm-meters
        mahalDist;  % Mahalanobis distance \Phi(G)
        likeProb;   % likelihood
        Cdi;        % inverse of data+forward modeling covariance matrix.
        %%%%%%%%%% The following are bounds on values %%%%%%%%%%
        lDepthMin;   % Minimum log-depth for a layer interface (besides top)
        lDepthMax;   % Maximum log-depth for a layer interface
        lRhoMin;     % Minimum log-resistivity
        lRhoMax;     % Maximum log-resistivity
        lHMin;       % Minimum log-thickness
        lVarMin;     % in LOG
        lVarMax;     % in LOG
        numMeasurements;
        %%%%%%%%%%%% How much to change values by %%%%%%%%%%%%%%
        lDepthChange;% standard deviation for depth changes in log space
        lRhoChange;  % standard deviation for rho changes in log space
        lVarChange;  % log, standard deviation for var changes
        badRunsThreshold; % how many bad runs you can have before it errors
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%% Core Functions %%%%%%%%%%
        
        %Constructor - initializes with a 1-layer model
        function obj = genericSln(pBounds,numMeasurements)
            %pBounds is a structure containing parameter bounds like min
            %and max values for depths, resistivity, variance; NOT in
            %log-space. numMeasurements is just the number of measurements
            %taken in the data set this is being used to invert.
            obj.numMeasurements = numMeasurements;
            obj.lDepths = [log10(0); nan*zeros(pBounds.maxLayers-1,1)];
            %depths,rhos must always be length numLayers, empty layers NaN
            obj.numLayers = nnz(~isnan(obj.lDepths));
            obj.var = pBounds.intlVar;
            obj.lDepthMin = log10(pBounds.depthMin);
            obj.lDepthMax = log10(pBounds.depthMax);
            obj.lRhoMin = log10(pBounds.rhoMin);
            obj.lRhoMax = log10(pBounds.rhoMax);
            obj.lRhos = [0.5*(obj.lRhoMin+obj.lRhoMax);... %average of bounds
                nan*zeros(pBounds.maxLayers-1,1)];
            obj.lVarMin = log10(pBounds.varMin);
            obj.lVarMax = log10(pBounds.varMax);
            obj.lHMin = (obj.lDepthMax - obj.lDepthMin)/...
                (2*pBounds.maxLayers); %Malinverno 2002, Appendix A1
            obj.lDepthChange = obj.lHMin;    % Update?
            obj.lRhoChange = obj.lRhoMax/2;  % Update?
            obj.lVarChange = log10(pBounds.varChange);
            obj.mahalDist = 0;
            obj.likeProb = 0;
            obj.misfit = 0;
            obj.Cdi = 0;
            obj.badRunsThreshold = ceil(log10(pBounds.numSteps)*20);
            %arbitrary
        end
        
        %This is meant to send the solution to another genericSln object
        function output = sendSln(obj)
            output.lDepths = obj.lDepths;
            output.lRhos = obj.lRhos;
            output.var = obj.var;
            output.misfit = obj.misfit;
            output.Cdi = obj.Cdi;
            output.residual = obj.residual;
        end
        
        %This is meant to recieve properties from another genericSln object
        function recieveSln(obj, input)
            obj.lDepths = input.lDepths;
            obj.lRhos = input.lRhos;
            obj.numLayers = nnz(~isnan(obj.lDepths));
            obj.Cdi = input.Cdi;
            obj.residual = input.residual;
            obj.var = input.var;
            obj.misfit = input.misfit;
            obj.calculatePosterior();
            if nnz(~isnan(obj.lDepths)) ~= nnz(~isnan(obj.lRhos))
                fprintf('Error: accepted solution # of layers dont match');
            end
        end
        
        %%%%%%%%%%%%%% Bound Checks %%%%%%%%%%%%%%%%%%%%
        
        % Sorts layers of proposed change. Should be called before
        % accepting sln or checking properties
        function sortLayers(obj)
            [obj.lDepths,idx]=sort(obj.lDepths);
            obj.lRhos = obj.lRhos(idx);
        end
        
        %Checks new depth values to make sure they're in bounds
        function good = checkDepthProperties(obj,proposedLDepth,indx)
            testLDepths = obj.lDepths(~isnan(obj.lDepths));
            %in addLayer, indx is numLayers+1, so indx<=length(testDepths)
            %only happens if coming from perturbDepth
            if indx <= length(testLDepths) %if perturbDepth
                testLDepths(indx) = []; %remove old depth before comparing
            end
            if (proposedLDepth > obj.lDepthMax-obj.lHMin ||...
                    proposedLDepth < obj.lDepthMin) || ...
                    (min(abs(proposedLDepth - testLDepths)) < obj.lHMin)
                %First set of criteria is checking depth bounds, second 
                %checks to make sure no layer is made too thin
                good = false;
            else
                good = true;
            end
        end
        
        %Checks to make sure rho values are in bounds
        function good = checkRhoProperties(obj,proposedLRho)
            if (proposedLRho < obj.lRhoMin || proposedLRho>obj.lRhoMax)
                good = false;
            else
                good = true;
            end
        end
        
        %Checks var values to be in bounds, assumes input in logspace
        function good = checkVarProperties(obj,input)
            if (input > obj.lVarMax || input < obj.lVarMin)
                good = false;
            else
                good = true;
            end
        end
        
        %%%%%%%%%% RANDOM WALK OPTIONS %%%%%%%%%%%%%%%
    %Note: once an option is chosen, it will attempt to complete that
    %option until it succeeds or it hits the number of attempts threshold
    %(badRunsThreshold). This is by design.

        %Alters proposed sln by changing layer depth, assumes >1 layer
        function perturbDepth(obj)
            success = false;
            nbad = 0;
            while ~success
                nbad=nbad+1;
                if(nbad>obj.badRunsThreshold)
                    error('nbad exceeded max,perturbDepth');
                end
                indx = randi([2,obj.numLayers]); %Don't touch first layer
                %Once layer is chosen, change depth and check to make sure
                %bounds aren't violated
                dummy = obj.lDepths(indx) + obj.lDepthChange*randn;
                success = obj.checkDepthProperties(dummy,indx);
            end
            obj.lDepths(indx) = dummy; %Make the change
            if ~issorted(obj.lDepths) %Layers may now be out of order
                obj.sortLayers();
            end
        end
        
        % Delete a layer, NOT the top layer
        %This concerns me and I should think more on it.
        function deleteLayer(obj)
            indx = randi([2,obj.numLayers]);
            obj.lDepths(indx:end) = [obj.lDepths(indx+1:end); NaN];%shift cells up
            obj.lRhos(indx:end) = [obj.lRhos(indx+1:end); NaN];
            obj.numLayers = nnz(~isnan(obj.lDepths));
        end
        
        %
        function addLayer(obj)
            success = false;
            nbad = 0;
            indx = obj.numLayers+1; %Index of the new layer
            %This does NOT mean newest layer will be lowest. depth is
            %randomly chosen, layers then sorted afterwards
            while ~success
                nbad=nbad+1;
                if(nbad>obj.badRunsThreshold)
                    error('nbad exceeded max,addLayer');
                end
                %Propose new depth and resistivity within bounds
                dummyLDepth = obj.lDepthMin +...
                    rand*(obj.lDepthMax - obj.lDepthMin);
                dummyLRho = obj.lRhoMin + rand*(obj.lRhoMax - obj.lRhoMin);
                if obj.checkDepthProperties(dummyLDepth,indx)
                    success = true;
                end
            end
            obj.lDepths(indx) = dummyLDepth;
            obj.lRhos(indx) = dummyLRho;
            obj.numLayers = nnz(~isnan(obj.lDepths));
            if ~issorted(obj.lDepths)
                obj.sortLayers();
            end
        end
        
        %
        function perturbRho(obj)
            success = false;
            nbad = 0;
            indx = randi([1,obj.numLayers]); %choose layer
            while ~success
                nbad=nbad+1;
                if(nbad>obj.badRunsThreshold)
                    error('nbad exceeded max ,perturbRho');
                end
                dummy = obj.lRhos(indx) + obj.lRhoChange*randn; %propose new rho
                success = obj.checkRhoProperties(dummy); %check
            end
            obj.lRhos(indx) = dummy; %actually change it
        end
        
        %
        function perturbVar(obj)
            %Var not in log-space so convert before updating
            success = false;
            nbad = 0;
            while ~success
                nbad = nbad+1;
                if (nbad>obj.badRunsThreshold)
                    error('nbad exceeded max, perturbVar')
                end
                lvar = log10(obj.var);
                dummy = lvar+(obj.lVarChange*randn);
                success = obj.checkVarProperties(dummy);
            end
            obj.var = 10^dummy;
        end
        
        %%%%%%%%%% GET IT %%%%%%%%%%%%%%%%%%%%
        %Extract properties for use outside of genericSln
        
        %Outputs depths and resistivities in not-log space
        function [depths,rhos] = getSolution(obj)
            depths = 10.^obj.lDepths;
            rhos = 10.^obj.lRhos;
        end
        
        %
        function N = getNumLayers(obj)
            N = obj.numLayers;
        end
        
        %
        function output = getMahalDist(obj)
            output = obj.mahalDist;
        end
        
        %
        function output = getVar(obj)
            output = obj.var; %remember var is NOT in log-space
        end
        
        function output = getMisfit(obj)
            output = obj.misfit;
        end
        
        function output = getLikeProb(obj)
            output = obj.likeProb;
        end
        
        function output = getCdi(obj)
            output = obj.Cdi;
        end
        
        
        %%%%%%%%%%% SET IT %%%%%%%%%%%%
        %for manually setting things from outside the object
        function setCdi(obj,inCdi)
            obj.Cdi = inCdi;
        end


        function setMisfit(obj,residual)
            % input should be (vector) residual and Cdi.
            obj.residual = residual;
            obj.misfit = norm(residual);
            obj.calculatePosterior();
        end
        
        function calculatePosterior(obj)
            obj.mahalDist = obj.residual'*obj.Cdi*obj.residual/obj.var;
            %Kolb_Lekic 2014 eq. 3 - ????
            obj.likeProb = exp(-0.5*obj.mahalDist)/...
                ((sqrt(2*pi*obj.var))^obj.numMeasurements);
            %Kolb Lekic 2014 eq 4
        end
    end
end