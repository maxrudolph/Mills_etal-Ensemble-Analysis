function data = createSyntheticData(measure, forwardModel)
%{
    Creates artificial model and synthetic measurements for
MillsSeniorThesisMain. The input structure, 'measure,' contains:
minDist, maxDist: scalars, min and max electrode spacing during 
measurements, in meters; numMeasurements: number of measurements taken;
noiseCoef: How noisy the measurements are; modelChoice: which model to use;
kMax: max number of layers (for inversion purposes). forwardModel is a
function handle for the choice of forward model being used (script)
    Creates: data a structure with fields x: values of electrode spacings
measurements; lambda: a 2D array based on x; fx: the model output with no
noise; y: model output with noise.
    Requires other scripts: modelGen, and something for forwardModel
    %}
    [trueDepths,trueRhos] = modelGen(measure.kMax,measure.modelChoice);
    minDist = log10(measure.minDist);
    maxDist = log10(measure.maxDist);
    data.x = logspace(minDist,maxDist,measure.numMeasurements); 
    %array of electrode spacings
    data.lambda = makeLambda(data.x);
    %lambda is an 11 x numMeasurements array of constants
    data.fx = forwardModel(trueDepths,trueRhos,data.lambda); 
    %The "true" model output given above parameters
    noiseVector = measure.noiseCoef.*data.fx.*randn(length(data.fx),1);
    %The absolute value of this shouldn't matter since it gets scaled by var.    
    data.y =  data.fx+noiseVector; %measurements with noise
    data.Cd = diag(data.y.^2);     % assume sig_f/f = const (constant relative error).
end
