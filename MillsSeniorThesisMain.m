clear
close all
rng(1) %seed to get same outcomes every time
%{
Inverse modeling 1D DC resistivity measurements. Based on tdmcmc_teaching
example. Use with calculateRho1D, genericMedium, and mcmcAlgorithm.
 Chris Mills 10/2020
%}

%% Set inversion options
options.kMax = 10; %max number of layers allowed in models
options.numSteps = 1e4; %total iterations for MCMC loop
options.saveStart = ceil(options.kMax*options.numSteps/...
    sum(2:options.kMax))/2;
% number of steps before end to start sampling. I have set it as is so that
% it does not start sampling until the max # of layers has been reached AND
% it has had time to test a few models with max # of layers
options.saveSkip = 10; %sample every (saveSkip)th step once sampling begins
options.samplePrior = false; % always accept proposed solution
options.intlVar = 1.0; %variance = how much misfit accepted.
options.alterVar = true; %If false, variance will never change

%% Define and measure starting model
%Measurement options
measure.minDist = 0.1; % Smallest electrode distance, meters
measure.maxDist = 1000; %  Largest electrode distance, meters
measure.numMeasurements = 21; %total # of measurements between minDist,maxDist
measure.noiseCoef = 0.1; %How "noisy" are the measurements, smaller values = less noisy
modelChoice = '4LayerA'; %choice of model

[trueDepths,trueRhos] = modelGen(options.kMax,modelChoice);

% Define true model
trueNumLayers = nnz(~isnan(trueDepths)); %k = number of layers
minDist = log10(measure.minDist);
maxDist = log10(measure.maxDist);
x = logspace(minDist,maxDist,measure.numMeasurements); 
measure.x = x;
%electrode spacings, meters
lambda = makeLambda(x); %lambda is an 11 x numMeasurements array of constants
measure.lambda = lambda;
fx = calculateRho1D(trueDepths,trueRhos,lambda); %The "true" model output given above parameters
noiseVector = measure.noiseCoef.*fx.*randn(length(fx),1);
% add normally distributed noise
y =  fx+noiseVector; %measurements with noise

%% Do the thing
ifSave = input('save?\n');

model = @(a,b,c) calculateRho1D(a,b,c);
results = mcmcAlgorithm( x,y,model,options ); %Do the inversion
if ifSave
    filename = ['Ensemble_', num2str(measure.noiseCoef), '_', choice,...
        '_', date, '.mat'];
    save(filename,'results','trueDepths','trueRhos','x','y','fx',...
        'options','noiseCoef');
end

%% Calculate and measure "mean" model
%Model represents the mean of all models (at each depth, the mean
%resistivity across all saved models for that depth)
%disp(['Creating mean model']);
%meanModel = createMeanModel(options.kMax,measure,results);

%% Plot Run Properties
disp(['Plotting']);
ensembleAnalysis(results,measure,model);

