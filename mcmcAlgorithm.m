function results = mcmcAlgorithm(data,model,options,pBounds)
%{
    Inputs:data.x and data.y are the observations (electrode spacings, 
apparent resistivities); model is handle to a function calculating 
forward-model given parameters, options contains fields kMax (max number of 
layers); numSteps (number of steps); mLPSCoefficient(int, involved with 
length of burn-in time); saveStart,saveSkip (ints controlling which runs to
save); samplePrior (bool, accept all proposed Slns or not); intlVar 
(variance); alterVar (bool, controls whether variance varies or not);
    Output: result is a structure with fields storedDepths; storedRhos;
storedSavedVars; storedLikelihoods; allMisfits; ensembleMisfits;
storedProbAccepts; numLayers; maxLayers; storedChoices;
    External scripts: genericMedium 
%}
%% Initializing
numMeasurements = length(data.x);
if options.alterVar %Whether or not variance can be altered
    randomOptions = 5;
else
    randomOptions = 4;
end

lambda = makeLambda(data.x);
totalSteps = options.numSteps; %total # of iterations to run
saveStart = totalSteps - options.saveStart; % what iteration to start saving
saveSkip = options.saveSkip;
pctSteps = round(totalSteps/100); %For printing status
numSavedRuns = ceil(options.saveStart/saveSkip); %How many runs to save

% Process covariance matrix
Cdi = pinv(data.Cd); % compute the Moore-Penrose pseudoinverse of the data covariance matrix

% Initialize solutions
layersProposed = genericMedium(pBounds,numMeasurements);
layersAccepted = genericMedium(pBounds,numMeasurements);

[depths,rhos] = layersAccepted.getSolution();
residual = data.y - model(depths,rhos,lambda);
layersAccepted.setCdi(Cdi);
layersProposed.setCdi(Cdi);
layersAccepted.setMisfit(residual);

% Pre-allocate memory for saving runs 
results.ensembleDepths = nan*zeros(pBounds.maxLayers,numSavedRuns);
results.ensembleRhos = nan*zeros(pBounds.maxLayers,numSavedRuns);
results.ensembleVars = zeros(1,numSavedRuns); %only for saved solutions
results.ensembleMisfits = zeros(1,numSavedRuns);
results.ensembleNumLayers = zeros(1,numSavedRuns);

results.allChoices = zeros(totalSteps,1);
results.allLikelihoods=zeros(totalSteps,1);
results.allMisfits = zeros(totalSteps,1);
results.allProbAccepts = zeros(totalSteps,1);
results.allVars = zeros(totalSteps,1);

% Slow down layer adding during the burn-in period
maxLayersPerStep = []; %This will be the max layers in any given step.
for iter = 2:pBounds.maxLayers
    maxLayersPerStep = [maxLayersPerStep;...
        iter*ones(options.mLPSCoefficient*iter,1)];
end
if length(maxLayersPerStep) < totalSteps
    maxLayersPerStep(end+1: totalSteps) = maxLayersPerStep(end);
end
results.maxLayers = maxLayersPerStep;
%% Main loop
saveStep = 0;
disp(['MCMC algorithm starting']);
for iter=1:totalSteps  %Number of steps in Markov Chain    
    resetProposedSln(layersAccepted,layersProposed);
    choice = chooseOption(layersProposed.getNumLayers(),...
        maxLayersPerStep(iter),randomOptions);
    switch choice
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
    end %end of switch statement
    
    if options.samplePrior
        probAccept = log(1);
    else
        [depths,rhos] = layersProposed.getSolution();
        residual = data.y - model(depths,rhos,lambda);       
        layersProposed.setMisfit(residual);
        %probAccept is calculated in ln space
        phi = layersAccepted.getMahalDist();
        phiPrime = layersProposed.getMahalDist();
        sigma2 = layersAccepted.getVar();
        sigma2Prime = layersProposed.getVar();
        k = layersAccepted.getNumLayers();
        kPrime = layersProposed.getNumLayers();
        probAccept = 0.5*(phi - phiPrime + numMeasurements*(log(sigma2) - ...
            log(sigma2Prime))) + log(k+1) - log(kPrime+1);
    end

    if ( isfinite(probAccept)&&( probAccept > log(rand)))
        acceptProposedSln(layersAccepted,layersProposed)
    end
    
    if round(iter/pctSteps)==iter/pctSteps
        %display status message
        disp(['MCMC algorithm ',num2str(100*(iter/totalSteps)),'% finished']);
    end
    
    if (iter>saveStart)&&(round(iter/saveSkip)==iter/saveSkip) %starts recording after ItsSaved models for every Inc_save-th model
        saveStep=saveStep+1;
        if saveStep == 1
            disp(['Saving begun']);
        end
        %Store properties for ensemble
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

%Choosing a random option
function output = chooseOption(currentLayers,maxLayers,numOptions)
    if currentLayers == 1 %If one layer model...
        % prob(death) = prob(perturbDepth) = 0;
        %numOptions = [4,5] => prob(addLayer) = [25%,20%]; 
        %prob(perturbRho) = [75,48%],prob(changeVar) = [0%,32%]
        tmp = rand(); %uniformly distributed
        if tmp<(1/numOptions)
            output = 3; %Add new layer.
        elseif (tmp < 1/numOptions + (1-1/numOptions)/(numOptions-3))
            output = 4; %perturb resistivity
        else
            output = 5; %var
        end
    elseif currentLayers >= maxLayers %If at max layers...
        %prob (addLayer) = 0; numOptions = [4,5] => 
        %prob(deleteLayer) = [25%,20%]; prob(perturbDepth) = [46.9%,37.3%];
        %prob(perturbRho) = [28.1%,31.3%], prob(changeVar) = [0%,11.4%]
        tmp = rand(); %uniformly distributed
        if tmp < 1/numOptions
            output = 2; %delete layer
        elseif (tmp < 1/numOptions + (1-1/numOptions)/(numOptions-2))
            output = 1; %perturb depth
        elseif (tmp < 1/numOptions + (2-2/numOptions)/(numOptions-2))
            output = 4; %perturb resistivity
        else
            output = 5;
        end
    else %normal circumstances
        output = randi([1,numOptions]);
    end
end
end