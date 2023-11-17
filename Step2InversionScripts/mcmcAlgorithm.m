function results = mcmcAlgorithm(data,model,options,pBounds)
%{
6/21/21
Performs an inversion by Markov-Chain Monte Carlo and generates a solution
ensemble. Inputs:
    data: structure with following fields (may have others, but only these
    are used):
        x: array of electrode spacings of virtual measurements
        y: measured apparent resistivities at each x, with noise added (the
            'actual data')
        lambda: lambda matrix associated with x, see makeLambda (needed for
            forward model)
        Cd: Covariance matrix
    model: a function handle for the forward model
    options: structure containing information governing the mcmc procedure
        numSteps: how many iterations of the main loop to take
        mLPSCoefficient: governs the length of the 'burn-in' period
        saveStart: How far from the end before sampling begins
        saveSkip: Once saving begins, save every saveSkipth solution
        alterVar: Controls whether hierarchical or not. 'True' = yes
            hierarchical and variance will be allowed to change. 'false' =
            no, constant variance maintained throughout.
        samplePrior: 'true' will sample prior distribution. Only set to
            true for testing purposes.
        pctSteps: Convenience, controls status update frequency
    pBounds: structure containing parameter bounds/information:
        intlVar: Initial variance
        maxLayers: maximum number of layers a sln can have.
        depthMin,depthMax,rhoMin,rhoMax,varMin,varMax: upper and lower
            bounds on depths and resistivities of layers in slns, and
            variance
        varChange: approximately controls how much the variance can shift
            if that is an option.
        numSteps: same as above, included here for the genericSln
            construction
        
    Output: results is a structure with fields:
        ensembleDepths,ensembleRhos,ensembleVars,ensembleMisfits,...
            ensembleNumLayers: Recording the properties of the ensemble
            members
        allChoices: records the choice picked at each step
        allLikelihoods,allMisfits,allProbAccepts,allVars,maxLayers: records
            these properties at each step (not just for ensemble members).
%}
%% Initializing
rng(1); %reproducibility
numMeasurements = length(data.x);
if options.alterVar %Sets whether or not variance can be altered
    randomOptions = 5;
else
    randomOptions = 4;
end

lambda = data.lambda;
totalSteps = options.numSteps; %total # of iterations to run
saveStart = totalSteps - options.saveStart; % what iteration to start saving
saveSkip = options.saveSkip;
pctSteps = round(totalSteps*options.pctSteps/100); %For printing status
numSavedRuns = ceil(options.saveStart/saveSkip); %How many runs to save

% Process covariance matrix
Cdi = pinv(data.Cd); % compute the Moore-Penrose pseudoinverse of the data
%covariance matrix

% Initialize 'solutions.' These are self-contained objects which have all
% the necessary protocols and information to perform loop steps.
%See 'genericSln' script
layersProposed = genericSln(pBounds,numMeasurements,Cdi);
layersAccepted = genericSln(pBounds,numMeasurements,Cdi);
[depths,rhos] = layersAccepted.getSolution();
acceptedGm = model(depths,rhos,lambda);
residual = data.y - acceptedGm;
layersAccepted.setMisfit(residual);
layersProposed.setMisfit(residual);
layersProposed.calculateRhoPrior();
layersAccepted.calculateRhoPrior();

% Pre-allocate memory for saving runs. Fields starting with 'ensemble...'
% only store information from saved solutions, which is toward the end of
% the loop. Fields starting with 'all' will catch from every step.
results.ensembleDepths = nan*zeros(pBounds.maxLayers,numSavedRuns);
results.ensembleRhos = nan*zeros(pBounds.maxLayers,numSavedRuns);
results.ensembleVars = zeros(1,numSavedRuns);
results.ensembleMisfits = zeros(1,numSavedRuns);
results.ensembleNumLayers = zeros(1,numSavedRuns);
results.ensembleGm = zeros(length(acceptedGm),numSavedRuns);

results.allChoices = zeros(totalSteps,1);
results.allLikelihoods=zeros(totalSteps,1);
results.allMisfits = zeros(totalSteps,1);
results.allProbAccepts = zeros(totalSteps,1);
results.allVars = zeros(totalSteps,1);

% Sets up a 'burn-in' period. Sets the maximum allowed # of layers for a
% sln at each step to gradually increase
maxLayersPerStep = []; %This will be the max layers in any given step.
for iter = 2:pBounds.maxLayers
    maxLayersPerStep = [maxLayersPerStep;...
        iter*ones(options.mLPSCoefficient*iter,1)];
end
if length(maxLayersPerStep) < totalSteps %this will most likely be true
    maxLayersPerStep(end+1:totalSteps) = maxLayersPerStep(end);
end
results.maxLayers = maxLayersPerStep;
%% Main loop
saveStep = 0;
disp(['MCMC algorithm starting']);
for iter=1:totalSteps  %Number of steps in Markov Chain
    resetProposedSln(layersAccepted,layersProposed);
    %Step 1: Set proposed sln = accepted sln
    choice = chooseOption(layersProposed.getNumLayers(),...
        maxLayersPerStep(iter),randomOptions);
    %Step 2: choose how proposed sln will be edited
    switch choice %Step 3: Edit proposed sln
        case 1 % Random option 1: Change the interface depth
            layersProposed.perturbDepth();
        case 2 %Random option 2: Delete a layer
            layersProposed.deleteLayer();
        case 3 % Random option 3: Add a new layer
            layersProposed.addLayer();
        case 4 %Random option 4: Change a layer's resistivity
            layersProposed.perturbRho();
        case 5 % Random option 5: Change noise variance
            layersProposed.perturbVar();
    end
    [depths,rhos] = layersProposed.getSolution();
    proposedGm = model(depths,rhos,lambda); % This is the forward model
    residual = data.y - proposedGm;
    layersProposed.setMisfit(residual);
    layersProposed.calculateRhoPrior();
    %Step 4: Run proposed sln through forward model to get misfit
    
    %Step 5: choose whether or not to accept the proposed sln
    %See my paper for this, or Malinverno 2002
    if options.samplePrior
        k = layersAccepted.getNumLayers();
        kPrime = layersProposed.getNumLayers();
        probAccept = log(k) - log(kPrime) + layersProposed.getPrior() - layersAccepted.getPrior();
    else
        %probAccept is calculated in ln space
        phi = layersAccepted.getMahalDist();
        phiPrime = layersProposed.getMahalDist();
        sigma2 = layersAccepted.getVar();
        sigma2Prime = layersProposed.getVar();
        k = layersAccepted.getNumLayers();
        kPrime = layersProposed.getNumLayers();
        probAccept = 0.5*(phi - phiPrime + numMeasurements*(log(sigma2) - ...
            log(sigma2Prime))) + log(k) - log(kPrime) + layersProposed.getPrior() - layersAccepted.getPrior();
        %Note that genericSln has a likeProb property, but it is preferable
        %to do this this way since its not the 'true' likeProb
    end
    
    if ( isfinite(probAccept)&&( probAccept > log(rand)))
        %Compare probAccept with a random number from uniform dist on (0,1)
        acceptProposedSln(layersAccepted,layersProposed);
        acceptedGm = proposedGm;
        %Proposed sln becomes new accepted sln
    end
    
    if round(iter/pctSteps)==iter/pctSteps
        %update user w/percent finished
        disp(['MCMC algorithm ',num2str(100*(iter/totalSteps)),'% finished']);
    end
    
    if (iter>saveStart)&&(round(iter/saveSkip)==iter/saveSkip)
        %starts saving slns once necessary criteria are met.
        saveStep=saveStep+1;
        if saveStep == 1
            disp(['Saving begun']);
        end
        %Store properties for ensemble
        results.ensembleG(:,saveStep) = acceptedGm;
        [results.ensembleDepths(:,saveStep),...
            results.ensembleRhos(:,saveStep)] = layersAccepted.getSolution();
        results.ensembleVars(saveStep) = layersAccepted.getVar();
        results.ensembleMisfits(saveStep) = layersAccepted.getMisfit();
    end
    %Store properties for entire run
    results.allMisfits(iter) = layersAccepted.getMisfit();
    results.allLikelihoods(iter) = layersAccepted.getLikeProb();
    results.allChoices(iter) = choice;
    results.allProbAccepts(iter) = probAccept;
    results.allVars(iter) = layersAccepted.getVar();
end

%% Wrap up
%save results to solution structure
results.ensembleNumLayers = sum(~isnan(results.ensembleDepths),1);

%% Functions
    function resetProposedSln(accepted,proposed)
        proposed.recieveSln(accepted.sendSln());
    end
%
    function acceptProposedSln(accepted,proposed)
        accepted.recieveSln(proposed.sendSln());
    end
end