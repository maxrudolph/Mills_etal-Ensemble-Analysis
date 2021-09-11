function filename = createSyntheticData(varargin)
%{
6/16/21
Generates synthetic data for subsequent inversion.
OPTIONAL inputs: noise coefficient (how noisy the data is, as a scalar) and
parameter subStructChoice, ie which subsurface structure to use (see
subStructGen script). Neither is necessary, so examples of how to call the 
script are: createSyntheticData, createSyntheticData(0.11),
createSyntheticData('subStructChoice','4LayerA'), or
createSyntheticData(0.11,'subStructChoice','4LayerA').
Saves in a file these things:
    forwardModel: Function handle of the forward model
    data: a structure containing data and related information
        x: values of electrode spacings, size numMeasurements, log-spaced
            from minDist to maxDist, in meters.
        lambda: a 2D array based on x, see makeLambda and calculateRho1D
        fx: the measurements with NO noise
        y: the measurements WITH noise
        Cd: covariance matrix, which in our case is
Then outputs a filename which can easily be fed into the inversion script,
which is the next step of the process.
%}
%% 0 Preliminary stuff ignore this section
addpath(genpath(fileparts(mfilename('fullpath'))))
%adds subfolders so you can use the scripts in them
defaultNoise = 0.11;
defaultSubStruct = '3LayerA';

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x>=0);
addOptional(p,'noise',defaultNoise,validScalarPosNum);
addParameter(p,'subStructChoice',defaultSubStruct,@ischar);
parse(p,varargin{:});


%% 1 User Set Options Here:

filterSize = 7;
%Choice of filter size for the forward model. Choices are 7,11,19. Based on
%Guptasarma 1982. This effects both which script is used for the forward
%model as well as the size of the lambda matrices, so it is an accuracy vs.
%computational expense tradeoff.

data.subStructChoice = p.Results.subStructChoice; 
%see subStructGen for choices

minDist = 0.1; maxDist = 1000; %in meters
numMeasurements = 21; %how many measurements
data.x = logspace(log10(minDist),log10(maxDist),numMeasurements);
% array of electrode spacings

data.noiseCoef = p.Results.noise; %How "noisy" are the measurements.

%% 2 Calculated Stuff:

switch filterSize
    case 7
        forwardModel = @(a,b,c) calculateRho1D07(a,b,c);
    case 11
        forwardModel = @(a,b,c) calculateRho1D11(a,b,c);
    case 19
        forwardModel = @(a,b,c) calculateRho1D19(a,b,c);
end %anything else is an invalid choice

[trueDepths,trueRhos] = subStructGen(data.subStructChoice);
data.lambda = makeLambda(data.x,filterSize); %lambda matrix for calculateRho1D
data.fx = forwardModel(trueDepths,trueRhos,data.lambda); %'true' output
rng(1); %re-seed random number generator for consistent noise pattern.
noiseVector = data.noiseCoef.*data.fx.*randn(length(data.fx),1);
%noise is added, Gaussian with mean 0 and std dev = noiseCoef*f(x)
data.y =  data.fx+noiseVector; %measurements with noise
data.Cd = diag(data.y.^2);
% assume sig_f/f = const (constant relative error).

%% 3 Save
filename = ['Data_', data.subStructChoice, '_',...
    num2str(data.noiseCoef), '_', date, '.mat'];
save(filename,'data','forwardModel');
end
