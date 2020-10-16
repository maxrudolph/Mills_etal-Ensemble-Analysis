function data = createSyntheticData(measure)
    [trueDepths,trueRhos] = modelGen(measure.kMax,measure.modelChoice);
    trueNumLayers = nnz(~isnan(trueDepths)); %k = number of layers
    minDist = log10(measure.minDist);
    maxDist = log10(measure.maxDist);
    data.x = logspace(minDist,maxDist,measure.numMeasurements); 
    %array of electrode spacings
    data.lambda = makeLambda(data.x);
    %lambda is an 11 x numMeasurements array of constants
    data.fx = calculateRho1D(trueDepths,trueRhos,data.lambda); 
    %The "true" model output given above parameters
    noiseVector = measure.noiseCoef.*data.fx.*randn(length(data.fx),1);
    % add normally distributed noise
    data.y =  data.fx+noiseVector; %measurements with noise
end
