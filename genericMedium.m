classdef genericMedium < handle
    %{
    Chris Mills 2021
    For use with the MCMC algorithm script. This object contains all the
    information for a proposed multi-layer medium containing some number of
    layers, each with a thickness (rather, a depth to the top of the layer)
    and a resistivity (rho). It also contains all information for updating
    itself and checking itself within bounds when a step is taken in the
    markov chain.
    %}
    properties (Access = private) %Fully self-contained
        %%%%%%%% Solution physical properties %%%%%%%%%%%%%%%%%%%%%%%%
        depths;     % current accepted solution log-depths, maxLayersx1
        rhos;       % current accepted solution log-resistivities,maxLayersx1
        numLayers;  % number of layers in proposed solution
        %%%%%%%%%%%%%%%%%% Solution math properties %%%%%%%%%%%%%%%%%%%%
        var;        % variance \sigma^2, NOT LOG
        residual;   % residual (vector)
        misfit;     % norm of residual
        mahalDist;  % Mahalanobis distance \Phi(G)
        likeProb;   %likelihood
        Cdi;        % inverse of data+forward modeling covariance matrix.
        %%%%%%%%%% The following are bounds on values %%%%%%%%%%
        depthMin;   % Minimum log-depth for a layer interface (besides top)
        depthMax;   % Maximum log-depth for a layer interface
        rhoMin;     % Minimum log-resistivity
        rhoMax;     % Maximum log-resistivity
        hMin;       % Minimum log-thickness
        varMin;     % in LOG
        varMax;     % in LOG
        numMeasurements;
        %%%%%%%%%%%% How much to change values by %%%%%%%%%%%%%%
        depthChange;% standard deviation for depth changes in log space
        rhoChange;  % standard deviation for rho changes in log space
        varChange;  % log, standard deviation for var changes 
        badRunsThreshold; % how many bad runs you can have before it errors
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods

%%%%%%%%%% Core Functions %%%%%%%%%%

        %Constructor - initializes with a 1-layer model
        function obj = genericMedium(pBounds,numMeasurements)
            %pBounds is a structure containing parameter bounds like min
            %and max values for depths, resistivity, variance; NOT in
            %log-space. numMeasurements is just the number of measurements
            %taken in the data set this is being used to invert.
            obj.numMeasurements = numMeasurements;
            obj.depths = [log10(0); nan*zeros(pBounds.maxLayers-1,1)];
            %depths,rhos must always be length numLayers, empty layers NaN
            obj.numLayers = nnz(~isnan(obj.depths));
            obj.setVar(pBounds.intlVar);
            obj.depthMin = log10(pBounds.depthMin);
            obj.depthMax = log10(pBounds.depthMax);
            obj.rhoMin = log10(pBounds.rhoMin);
            obj.rhoMax = log10(pBounds.rhoMax);
            obj.rhos = [obj.rhoMin+(0.5*(obj.rhoMax - obj.rhoMin));...
                nan*zeros(pBounds.maxLayers-1,1)];
            obj.varMin = log10(pBounds.varMin);
            obj.varMax = log10(pBounds.varMax);
            obj.hMin = (obj.depthMax - obj.depthMin)/(2*pBounds.maxLayers);
            %Malinverno 2002, Appendix A1
            obj.depthChange = obj.hMin;    % Update?
            obj.rhoChange = obj.rhoMax/2;  % Update?
            obj.varChange = log10(pBounds.varChange);
            obj.mahalDist = 0;
            obj.likeProb = 0;
            obj.misfit = 0;
            obj.badRunsThreshold = ceil(log10(pBounds.numSteps)*20);
            %arbitrary
        end
        
        %This is meant to send the solution to another genericMedium object
        function output = sendSln(obj)
            output.depths = obj.depths;
            output.rhos = obj.rhos;
            output.var = obj.var;
            output.misfit = obj.misfit;
            output.Cdi = obj.Cdi;
            output.residual = obj.residual;
        end
        
        %This is meant to recieve properties from another genericMedium
        %object
        function recieveSln(obj, input)
            obj.depths = input.depths;
            obj.rhos = input.rhos;
            obj.numLayers = nnz(~isnan(obj.depths));
            obj.Cdi = input.Cdi;
            obj.residual = input.residual;
            obj.setVar(input.var);
            obj.misfit = input.misfit;
            obj.calculatePosterior();
            if nnz(~isnan(obj.depths)) ~= nnz(~isnan(obj.rhos))
                fprintf('Error: accepted solution # of layers dont match');
            end
        end
        
%%%%%%%%%%%%%% Checks %%%%%%%%%%%%%%%%%%%%        
        
        %
        function sortLayers(obj)
            %sorts layers of proposed solution. Should be called before
            %accepting solution or checking properties
            [~,idx]=sort(obj.depths);
            obj.rhos = obj.rhos(idx);
            obj.depths = obj.depths(idx);
        end
        
        %
        function good = checkDepthProperties(obj,proposedDepth,indx)
            %checks new depth values to make sure in bounds
            testDepths = obj.depths(~isnan(obj.depths));
            %in addLayer, indx is numLayers+1, so indx<=length(testDepths)
            %only happens if coming from perturbDepth
            if indx <= length(testDepths) %if perturbDepth
                testDepths(indx) = []; %remove old layer depth before comparison
            end
            if (proposedDepth > obj.depthMax-obj.hMin || proposedDepth < obj.depthMin) || ...
                    (min(abs(proposedDepth - testDepths)) < obj.hMin)
                %First is checking depth limits, second checks to make sure
                %no layer is made too thin
                good = false;
            else
                good = true;
            end
        end
        
        %
        function good = checkRhoProperties(obj,proposedRho)
            %checks to make sure rho values are in bounds
            if (proposedRho < obj.rhoMin || proposedRho>obj.rhoMax)
                good = false;
            else
                good = true;
            end
        end
        
        %
        function good = checkVarProperties(obj,input)
            %checks var values to be in bounds
            %assumes input in logspace
            if (input > obj.varMax || input < obj.varMin)
                good = false;
            else
                good = true;
            end
        end
        
%%%%%%%%%% RANDOM WALK OPTIONS %%%%%%%%%%%%%%%
        
        %
        function perturbDepth(obj)
            %Alters proposed solution by changing a layer depth
            %Assumes more than a one-layer model (don't use otherwise)
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
                dummy = obj.depths(indx) + obj.depthChange*randn; 
                success = obj.checkDepthProperties(dummy,indx);
            end
            obj.depths(indx) = dummy; %Make the change
            if ~issorted(obj.depths)
                obj.sortLayers(); %Layers may be out of order now
            end
        end
        
        %
        function deleteLayer(obj)
            %The top layer can't be deleted
            indx = randi([2,obj.numLayers]);
            obj.depths(indx:end) = [obj.depths(indx+1:end); NaN];%shift cells up
            obj.rhos(indx:end) = [obj.rhos(indx+1:end); NaN];
            obj.numLayers = nnz(~isnan(obj.depths));
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
                dummyDepth = obj.depthMin + rand*(obj.depthMax - obj.depthMin);
                dummyRho = obj.rhoMin + rand*(obj.rhoMax - obj.rhoMin);
                %Next step really isn't necessary I think but I'm paranoid
                if obj.checkDepthProperties(dummyDepth,indx) &&...
                        obj.checkRhoProperties(dummyRho)
                    success = true;
                end
            end
            obj.depths(indx) = dummyDepth;
            obj.rhos(indx) = dummyRho;
            obj.numLayers = nnz(~isnan(obj.depths)); %replace with +1?
            if ~issorted(obj.depths)
                obj.sortLayers();
            end
        end
        
        %
        function perturbRho(obj)
            success = false;
            nbad = 0;
            while ~success
                nbad=nbad+1;
                if(nbad>obj.badRunsThreshold)
                    error('nbad exceeded max ,perturbRho');
                end
                indx = randi([1,obj.numLayers]); %choose layer
                dummy = obj.rhos(indx) + obj.rhoChange*randn; %propose new rho
                success = obj.checkRhoProperties(dummy); %check
            end
            obj.rhos(indx) = dummy; %actually change it
        end
        
        %
        function perturbVar(obj)
        %Var not in log-space,but varChange is, so convert before updating
            success = false;
            nbad = 0;
            while ~success
                nbad = nbad+1;
                if (nbad>obj.badRunsThreshold)
                    error('nbad exceeded max, perturbVar')
                end
                dvar = log10(obj.var);
                dummy = dvar+(obj.varChange*randn);
                success = obj.checkVarProperties(dummy);
            end
            obj.setVar(10^dummy);
        end
          
%%%%%%%%%% GET IT %%%%%%%%%%%%%%%%%%%%
%all should be pretty self-explanatory
        
        %Outputs depths and resistivities in not-log space
        function [depths,rhos] = getSolution(obj)
            depths = 10.^obj.depths;
            rhos = 10.^obj.rhos;
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
        
        
%%%%%%%%%%% SET IT %%%%%%%%%%%%%%%% %for manually setting stuff
        function setCdi(obj,Cdi)
            obj.Cdi = Cdi;
        end
        
        function setVar(obj,input)
            obj.var = input;
        end
        
        function setMisfit(obj,residual)
            % input should be (vector) residual and Cdi.
            obj.residual = residual;
            obj.misfit = norm(residual);
            obj.calculatePosterior();
        end
        
        function calculatePosterior(obj)
            obj.mahalDist = obj.residual'*obj.Cdi*obj.residual/obj.var;
            %Kolb_Lekic 2014 eq. 3
            obj.likeProb = exp(-0.5*obj.mahalDist)/...
                ((sqrt(2*pi*obj.var))^obj.numMeasurements);
            %Kolb Lekic 2014 eq 4
        end
    end
end