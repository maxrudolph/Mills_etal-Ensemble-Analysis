clear
close all
rng(1) %seed to get same outcomes every time
%{
Inverse modeling 1D DC resistivity measurements. Based on tdmcmc_teaching
example. This script is for doing an inverion to get an ensemble solution,
use ensembleAnalysis to analyze an already existing ensemble.
Requires other scripts: createSyntheticData, mcmcAlgorithm, 
ensembleAnalysis, calculateRho1D.
Chris Mills 10/2020
%}

forwardModel = @(a,b,c) calculateRho1D(a,b,c);

%% Step 1: Load data or create model/generate synthetic measurements
%{
ifLoadData = false;%input('Load data? true/false\n');

if ifLoadData
    % ...Not setup yet to take real data...
else
    %Artifical model setup:
    measure.modelChoice = '3LayerA'; %currently setup: 3LayerA, 4LayerA
    %Measurement options
    measure.minDist = 0.1; % Smallest electrode distance, meters
    measure.maxDist = 1000; %  Largest electrode distance, meters
    measure.numMeasurements = 21; %total # of measurements
    measure.noiseCoef = 0.1; %How "noisy" are the measurements
end
%}
%% Set inversion options
options.kMax = 10; %max number of layers allowed in models
options.numSteps = 1e4; %total iterations for MCMC loop. 1e7+ recommended
options.mLPSCoefficient = 1e5;
%mLPS = max layers per step. Set higher for longer 'burn-in' period.
options.saveStart = floor(options.numSteps/2);
%saveStart is the # of steps before end to start sampling. Should not
%sample until max # of layers has been reached AND it has had time to test
%several models with max # of layers.
options.saveSkip = 10; %sample every (saveSkip)th step once sampling begins
options.samplePrior = false; % true = always accept proposed solution
options.intlVar = 1.0; %variance = how much misfit accepted.
options.alterVar = true; %If false, model variance will never change

%If numSteps is not sufficiently high, reset above quantities appropriately
if options.mLPSCoefficient*sum(2:options.kMax) > options.numSteps
    options.mLPSCoefficient = floor(options.numSteps/sum(2:options.kMax));
    disp(['Not enough steps, changing burn-in time']);
end
if options.numSteps - options.saveStart < options.mLPSCoefficient*...
        (sum(2:options.kMax-1)+(0.5*options.kMax))
    options.saveStart = options.numSteps - (options.mLPSCoefficient*...
        ceil(sum(2:options.kMax-1)+(0.5*options.kMax)));
    disp(['Changing saveStart time']);
end

data.x = feetToMeters([10,25,40,50,100,125,150]);
data.y = [7177.30,4230.1,6228.1,8701.7,591356.0,8036.3,...
    99022.1];
data.lambda = makeLambda(data.x);
%data = createSyntheticData(measure, forwardModel); %creates measurements

%% Do inversion
results = mcmcAlgorithm(data,forwardModel,options); %Do the inversion

%{
noiseCoefValues = [0.01,0.05,0.1,0.15,0.2];

parfor irun = 1:length(noiseCoefValues)
    thisMeasure = measure;
    thisMeasure.noiseCoef = noiseCoefValues(irun);
    data = createSyntheticData(thisMeasure, forwardModel); %creates measurements
    results = mcmcAlgorithm(data,forwardModel,options);
    filename = ['Ensemble_', thisMeasure.modelChoice, '_',...
        num2str(thisMeasure.noiseCoef), '_', date, '.mat'];
    doSaving(filename,results,data,thisMeasure,options,forwardModel);
end
%}

%% Plot Run Properties
disp(['Plotting']);


ifSave = true;%false;%input('save? true/false\n');
if ifSave
    filename = ['Ensemble_ActualData_',date, '.mat'];
    save(filename,'results','data','options','forwardModel');
    ensembleAnalysis2(filename);
end
%}