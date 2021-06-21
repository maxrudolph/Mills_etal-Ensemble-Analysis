function filename = createSyntheticData
%{
6/16/21
Generates synthetic data for subsequent inversion.       
Outputs:
    forwardModel: Function handle of the forward model
    measure: a structure containing measurement details - 
        minDist, maxDist: scalars representing minimum and maximum
            electrode spacings, in meters
        numMeasurements: integer number of electrode spacings
        noiseCoef: How noisy you want the data to be, in the form of a
            coefficient. 0 is no noise, 0.2 is very noisy, etc.
        subStructChoice: string defining the subsurface structure. see
            subStructGen, if you want custom structure then add it there
    data: a structure containing data and related information
        x: values of electrode spacings, size numMeasurements, log-spaced 
            from minDist to maxDist, in meters.
        lambda: a 2D array based on x, see makeLambda and calculateRho1D
        fx: the measurements with NO noise
        y: the measurements WITH noise
        Cd: covariance matrix
%}
addpath(genpath(fileparts(mfilename('fullpath'))))
%adds subfolders so you can use the scripts in them


%% User Set Options:

forwardModel = @(a,b,c) calculateRho1D(a,b,c);

data.subStructChoice = '3LayerA'; %see subStructGen for choices

minDist = 0.1; maxDist = 1000; %in meters
numMeasurements = 21; %how many measurements
data.x = logspace(log10(minDist),log10(maxDist),numMeasurements);
% array of electrode spacings

data.noiseCoef = 0.11; %How "noisy" are the measurements.

%% Calculated Stuff:

[trueDepths,trueRhos] = subStructGen(data.subStructChoice);
data.lambda = makeLambda(data.x); %lambda matrix for calculateRho1D
data.fx = forwardModel(trueDepths,trueRhos,data.lambda); %'true' output
noiseVector = data.noiseCoef.*data.fx.*randn(length(data.fx),1);
%noise is added, Gaussian with mean 0 and std dev = noiseCoef*f(x)
data.y =  data.fx+noiseVector; %measurements with noise
data.Cd = diag(data.y.^2);
% assume sig_f/f = const (constant relative error).

%% Save
filename = ['Data_', data.subStructChoice, '_',...
    num2str(data.noiseCoef), '_', date, '.mat'];
save(filename,'data','forwardModel');
end
