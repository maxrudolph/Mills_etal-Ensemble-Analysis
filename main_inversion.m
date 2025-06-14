clear;
close all;
delete(gcp('nocreate'));
parpool(3);

%function saveEnsemblesLoop
noiseLevels = [0,0.01,0.02,0.05,0.1,0.2];
% subStructs = {'1LayerA','3LayerA','4LayerA'};
subStructs = {'3LayerA'};
nstruct = length(subStructs);
nnoise = length(noiseLevels);
noiseLevels = repmat(noiseLevels,1,nstruct)';
subStructs = repmat(subStructs,nnoise,1);
% subStructs = vertcat(subStructs{:,1},subStructs{:,2},subStructs{:,3});
subStructs = {subStructs{:}}';

mex ./Step1DataGenerationScripts/calculateRho1D11_mex.c

%Create data files
tic
parfor i = 1:length(noiseLevels)
    if i>length(noiseLevels) %sample the prior for one subStruct
        a = createSyntheticData(0,'subStructChoice','3LayerA');
        b = inversion(a,'priorOn',true);
    else
        a = createSyntheticData(noiseLevels(i),'subStructChoice',subStructs{i});
        b = inversion(a,'piecewiseLinear',true);
    end
end
toc

