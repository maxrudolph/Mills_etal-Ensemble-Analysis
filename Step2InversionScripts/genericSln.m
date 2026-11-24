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
    
%"l-" means log, so "lDepths" is in log-space, etc. Keeps track of what's 
%log and what isn't
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
        prior;      % prior probability for the model
    end
    
    properties(SetAccess = immutable) %Should not be changed once defined
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
        %% The type of prior used.
        priorChoice;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%% Core Functions %%%%%%%%%%
        
        %Constructor - initializes with a 1-layer model
        function obj = genericSln(pBounds,numMeasurements,Cdi)
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
            obj.lHMin = log10(pBounds.hMin);
            obj.lDepthChange = log10(pBounds.depthChange);
            obj.lRhoChange = log10(pBounds.rhoChange);         
            obj.lVarChange = log10(pBounds.varChange);
            obj.mahalDist = 0;
            obj.likeProb = 0;
            obj.misfit = 0;
            obj.Cdi = Cdi;
            obj.badRunsThreshold = ceil(log10(pBounds.numSteps)*20);
            obj.priorChoice = pBounds.priorChoice;
            %arbitrary
        end
        
        %This is meant to send the solution to another genericSln object
        function output = sendSln(obj)
            output.lDepths = obj.lDepths;
            output.lRhos = obj.lRhos;
            output.var = obj.var;
            output.misfit = obj.misfit;
            output.residual = obj.residual;
        end
        
        %This is meant to recieve properties from another genericSln object
        function recieveSln(obj, input)
            obj.lDepths = input.lDepths;
            obj.lRhos = input.lRhos;
            obj.numLayers = nnz(~isnan(obj.lDepths));
            obj.residual = input.residual;
            obj.var = input.var;
            obj.misfit = input.misfit;
            obj.calculatePosterior();
            if nnz(~isnan(obj.lDepths)) ~= nnz(~isnan(obj.lRhos))
                fprintf('Error: accepted solution # of layers dont match');
            end
        end
        
        %%%%%%%%%%%%%% Bound Checks %%%%%%%%%%%%%%%%%%%%
        
        % Compute available depth intervals
        function [validDepths,availDepth] = getValidDepths(obj)
            % This function computes the available depth ranges.
            % validDepths is a 2xN vector that contains the starting and
            % ending depth values of each available range.
            % availDepth is the total log10-depth range in which a layer
            % interface may be placed.
            %
            % if there are k layers, there are k-1 control points in the
            % valid depth range and therefore up to (k-1)+1 valid ranges
            % for placement
            if obj.numLayers == 1
                validDepths = {[obj.lDepthMin obj.lDepthMax]};
                availDepth = obj.lDepthMax-obj.lDepthMin;
            else                
                validDepths = zeros(2,obj.numLayers);
                ind=1;
                if obj.lDepths(2)-obj.lHMin < obj.lDepthMin
                    % range overlaps with start of interval  
                    
                else
                    % range does not overlap with start of interval
                    validDepths(:,ind) = [obj.lDepthMin, obj.lDepths(2)-obj.lHMin]; ind=ind+1;
                end
                for i=3:(obj.numLayers-1)
                    validDepths(:,ind) = [obj.lDepths(i-1)+obj.lHMin, obj.lDepths(i)-obj.lHMin]; ind=ind+1;
                end
                % for last layer - consider whether it overlaps with
                % endpoint of interval.
                i=obj.numLayers;
                if obj.lDepths(i)+obj.lHMin > obj.lDepthMax
                    % no interval below this depth
                else
                    validDepths(:,ind) = [obj.lDepths(i)+obj.lHMin, obj.lDepthMax]; ind = ind+1;
                end
                validDepths = validDepths(:,1:ind-1);
                availDepth = sum( validDepths(2,:)-validDepths(1,:) );
            end          
        end

        % function perturbWithinValidRange(obj,zstart,dz,validDepths)
        %     % this function computes a perturbation of depth starting at
        %     % zstart and ending at dz, skipping over any excluded zones
        %     % due to the presence of other control points.
        %     % find starting 'cell' within valid depths
        %     ind = find(zstart > validDepths(2,:),1,'first');
        %     distance_remaining = dz;
        %     while 
        % 
        % 
        % 
        % end

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
            if (proposedLDepth > obj.lDepthMax ||...
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
        
        % Compute the prior on rho
        % Case 1 - assumes a flat prior on Rho
        % Case 2 - assumes a gaussian prior centered on log10(rho) = 3
        % with sigma = 1 (in log rho units)
        function calculatePrior(obj)                 
            % switch obj.rhoPrior
            % case 1
            % obj.prior = 0;
            % case 2              
            
            % prior on each layer rho is normal distribution with
            % sigma=1 centered at log10(rho) = 3                                       
            % achieved close to the right distribution:
            if ~obj.checkVarProperties(log10(obj.var))
                obj.prior = log(0);
            elseif min(obj.lRhos(1:obj.numLayers)) < obj.lRhoMin || max(obj.lRhos(1:obj.numLayers)) > obj.lRhoMax
                obj.prior = log(0);
            elseif min(diff(obj.lDepths(1:obj.numLayers))) < obj.lHMin
                obj.prior = log(0);
            else
                switch obj.priorChoice
                    case 1                     
                        obj.prior = -log(obj.numLayers); %prior on k is 1/k
                    case 2
                        % flat prior on number of layers. Prior on rho gets
                        % wrapped up in proposal ratios below.
                        obj.prior = -log(obj.numLayers);
                end
            end
        end

        %%%%%%%%%% RANDOM WALK OPTIONS %%%%%%%%%%%%%%%
        %Note: once an option is chosen, it will attempt to complete that
        %option until it succeeds or it hits the number of attempts threshold
        %(badRunsThreshold). This is by design.
        
        %Alters proposed sln by changing layer depth, assumes >1 layer
        function priorRatio = perturbDepth(obj)
            success = false;
            nbad = 0;

            % [oldValidRanges,oldValidDepth] = obj.getValidDepths();

            while ~success
                nbad=nbad+1;
                if(nbad>obj.badRunsThreshold)
                    error('nbad exceeded max,perturbDepth');
                end
                indx = randi([2,obj.numLayers]); %Don't touch first layer
                %Once layer is chosen, change depth and check to make sure
                %bounds aren't violated
                dummy = obj.lDepths(indx) + obj.lDepthChange*(2*rand()-1);
                if dummy > obj.lDepthMax || dummy < obj.lDepthMin
                    priorRatio = log(0.0); % make a proposal but with zero probability of acceptance
                    break;
                else
                    layer_indices = 2:obj.numLayers;
                    not_perturbed = setdiff(layer_indices,indx); % list of unperturbed layers.
                    % [validRanges,validDepth] = obj.getValidDepths();

                    % if isempty(not_perturbed)
                    %     priorRatio = log(1.0);
                    %     break
                    if min( abs(obj.lDepths(not_perturbed)-dummy)) <= obj.lHMin % proposed layer is too thin.
                        priorRatio = log(0.0);
                        break                    
                    else
                        priorRatio = log(1.0);
                        break
                    end
                end
            end
            % success = obj.checkDepthProperties(dummy,indx);
            %end
            obj.lDepths(indx) = dummy; %Make the change
            if ~issorted(obj.lDepths) %Layers may now be out of order
                obj.sortLayers();
            end
        end
        
        % Delete a layer, NOT the top layer
        function priorRatio = deleteLayer(obj)
            indx = randi([2,obj.numLayers]);
            rhotmp = obj.lRhos(indx); % log-resistivity of the killed layer
            obj.lDepths(indx:end) = [obj.lDepths(indx+1:end); NaN];%shift cells up
            obj.lRhos(indx:end) = [obj.lRhos(indx+1:end); NaN];
            % obj.lRhos(indx) = 0.5*(obj.lRhos(indx) + rhotmp); % average per Malinverno 2002 Appendix.
            obj.numLayers = nnz(~isnan(obj.lDepths));
            if indx == 1
                obj.lDepths(1) = log(0);
            end
            % success = true;
            if obj.priorChoice == 1
                priorRatio = log(1.0);
            else
                sig_rho = 1.0;
                rho_bar = 3.0;
                priorRatio = 0.5*log(2*pi*sig_rho^2) - log(obj.lRhoMax-obj.lRhoMin) + 0.5/sig_rho^2*(rhotmp-rho_bar)^2;
            end
        end
        
        %
        function priorRatio = addLayer(obj)
            success = false;
            nbad = 0;
            indx = obj.numLayers+1; %Index of the new layer
            %This does NOT mean newest layer will be lowest. depth is
            %randomly chosen, layers then sorted afterwards
            
            while ~success % keep trying depths until we succeeed in adding a layer
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
            if obj.priorChoice == 1
                priorRatio = log(1.0);
            else
                sig_rho = 1.0;
                rho_bar = 3.0;
                priorRatio = log(obj.lRhoMax-obj.lRhoMin) - 0.5*log(2*pi*sig_rho^2) - 0.5*(dummyLRho-rho_bar)^2/sig_rho^2;
            end

        end
        
        %
        function priorRatio = perturbRho(obj)
            % success = false;
            % nbad = 0;
            indx = randi([1,obj.numLayers]); %choose layer
            % while ~success
                % nbad=nbad+1;
                % if(nbad>obj.badRunsThreshold)
                    % error('nbad exceeded max ,perturbRho');
                % end
                % note - this assigns rho from a normal distribution with
                % 95% CI between 1-5
                % dummy = 3 + randn;% sigma = 1, 2*sigma=2...
                % this assigns rho from a flat prior.
                oldRho = obj.lRhos(indx);
                dummy = obj.lRhos(indx) + obj.lRhoChange*randn; %propose new rho
                success = obj.checkRhoProperties(dummy); %check
            % end
            obj.lRhos(indx) = dummy; %actually change it
            if success
                if obj.priorChoice == 1
                    priorRatio = log(1);
                else
                    rho_bar = 3.0;
                    sig_rho = 1.0;
                    priorRatio = -0.5*(dummy-rho_bar)^2/sig_rho^2 - -0.5*(oldRho-rho_bar)^2/sig_rho^2;
                end
            else
                priorRatio = log(0);
            end
        end
        
        %
        function priorRatio = perturbVar(obj)
            %Var not in log-space so convert before updating
            % success = false;
            % nbad = 0;
            % while ~success
                % nbad = nbad+1;
                % if (nbad>obj.badRunsThreshold)
                    % error('nbad exceeded max, perturbVar')
                % end
                lvar = log10(obj.var);
                dummy = lvar+(obj.lVarChange*randn);
                success = obj.checkVarProperties(dummy);
            % end
            obj.var = 10^dummy;
            if success
                priorRatio = log(1);
            else 
                priorRatio = log(0);
            end
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
        
        function output = getPrior(obj)
            obj.calculatePrior();
            output = obj.prior;
        end
        
        %%%%%%%%%%% SET IT %%%%%%%%%%%%
        %can be used to manually set things from outside the object      
        
        function setMisfit(obj,residual)
            % input should be (vector) residual and Cdi.
            obj.residual = residual;
            obj.misfit = norm(residual);
            obj.calculatePosterior();
        end
        
        function calculatePosterior(obj)
            obj.mahalDist = obj.residual'*obj.Cdi*obj.residual/obj.var;
            %Kolb_Lekic 2014 eq. 3, divided by variance, because of the way
            %we set up the covariance matrix 
            obj.likeProb = exp(-0.5*obj.mahalDist)/...
                ((sqrt(2*pi*obj.var))^obj.numMeasurements);
            %Kolb Lekic 2014 eq 4. Note that this isn't the actual
            %likelihood probability, we leave out the determinant of the
            %covariance matrix which is effectively a constant
            %Note that likeProb is not used in any calculations, it is just
            %calculated for storage. 
        end
    end
end