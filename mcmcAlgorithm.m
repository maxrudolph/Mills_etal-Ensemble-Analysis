function result = mcmcAlgorithm(data,model,options)
%{
This routine is meant to be used with MillsSeniorThesisMain, 
 calculateRho1D, and the genericMedium class definition.
 xObserved and yObserved are the observations (electrode spacings, apparent
resistivities); model is handle to a function calculating forward-model 
 given parameters, options contains fields specifying max number of 
layers, number of steps, steps to save, and whether to sample prior pdf
%}
%% Initializing
% Create bounds on parameter values. These bounds are based on Appendix A
% in Malinverno 2002. See also the "genericMedium" constructor function
%Bound parameters. Bounds based on Appendix A, Malinverno 2002
numMeasurements = length(data.x);
pBounds.maxLayers = options.kMax; % max # of layers in a given model
pBounds.depthMin = min(data.x); %min depth for layer interface (not top)
pBounds.depthMax = max(data.x); % max depth for layer interface
pBounds.rhoMin = 1e-1; % min resistivity, NEEDS UPDATE
pBounds.rhoMax = 1e8; % max resistivity, NEEDS UPDATE
pBounds.varMin = 1e-6; % valid?  
pBounds.varMax = 1e3; % valid?
pBounds.varChange = 1e-1;  %valid?
pBounds.intlVar = options.intlVar; %initial variance
if options.alterVar
    randomOptions = 5;
else
    randomOptions = 4;
end
abciss = [-0.420625;-0.20265625;0.0153125;0.23328125;0.45125;...
    0.66921875;0.8871875;1.10515625;1.323125;1.54109375;1.7590625];

[xGrid,abcissGrid] = meshgrid(data.x,abciss);
lambda = 10.^(abcissGrid-log10(xGrid));

totalSteps = options.numSteps; %total # of iterations to run
saveStart = totalSteps - options.saveStart; % what iteration to start saving
saveSkip = options.saveSkip;
pctSteps = round(totalSteps/100); %For printing status
numSavedRuns = ceil(options.saveStart/saveSkip); %How many runs to save

% Initialize solutions
layersProposed = genericMedium(pBounds,numMeasurements);
layersAccepted = genericMedium(pBounds,numMeasurements);

[depths,rhos] = layersAccepted.getSln();
layersAccepted.setMisfit(data.y - model(depths,rhos,lambda));

% Pre-allocate memory for saving runs 
result.storedDepths = nan*zeros(pBounds.maxLayers,numSavedRuns);
result.storedRhos = nan*zeros(pBounds.maxLayers,numSavedRuns);
result.storedSavedVars = zeros(1,numSavedRuns); %only for saved solutions
result.storedLikelihoods=zeros(totalSteps,1);
result.allMisfits = zeros(totalSteps,1);
result.ensembleMisfits = zeros(1,numSavedRuns);
%result.storedChoices = zeros(totalSteps,1);
result.storedProbAccepts = zeros(totalSteps,1);
result.numLayers = zeros(1,numSavedRuns);

% Slow down layer adding during the burn-in period
maxLayersPerStep = []; %This will be the max layers in any given step.
mLPSCoefficient = 1e4;
if mLPSCoefficient*sum(2:pBounds.maxLayers) > totalSteps
    mLPSCoefficient = ceil(totalSteps/sum(2:pBounds.maxLayers));
    disp(['Not enough steps, changing burn-in time']);
end
for iter = 2:pBounds.maxLayers
    maxLayersPerStep = [maxLayersPerStep; iter*ones(mLPSCoefficient*iter,1)];
end
maxLayersPerStep(end+1:totalSteps) = maxLayersPerStep(end);
result.maxLayers = maxLayersPerStep;

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
        [depths,rhos] = layersProposed.getSln();
        layersProposed.setMisfit(data.y - model(depths,rhos,lambda));
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
    
    if (iter>saveStart)&&(round(iter/saveSkip)==iter/saveSkip)  %starts recording after ItsSaved models for every Inc_save-th model
        saveStep=saveStep+1;
        if saveStep == 1
            disp(['Saving begun']);
        end
        %Store properties for ensemble
        [result.storedDepths(:,saveStep),...
            result.storedRhos(:,saveStep)] = layersAccepted.getSln();
        result.storedSavedVars(saveStep) = layersAccepted.getVar();
        result.ensembleMisfits(saveStep) = layersAccepted.getMisfit(); 
    end
    %Store properties for entire run
    result.allMisfits(iter) = layersAccepted.getMisfit(); 
    result.storedLikelihoods(iter) = layersAccepted.getLikeProb();
    result.storedChoices(iter) = choice;
    result.storedProbAccepts(iter) = probAccept;
end

%% Wrap up
%save results to solution structure
result.numLayers = sum(~isnan(result.storedDepths),1);
%% Functions
function resetProposedSln(accepted,proposed)
    proposed.recieveSln(accepted.sendSln());
end

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