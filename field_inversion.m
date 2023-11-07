clear;
close all;

%function saveEnsemblesLoop
noiseLevels = [0]; % no additional noise.
subStructs = {'field'};
nstruct = length(subStructs);
nnoise = length(noiseLevels);
noiseLevels = repmat(noiseLevels,1,nstruct)';
subStructs = repmat(subStructs,nnoise,1);
% subStructs = vertcat(subStructs{:,1},subStructs{:,2},subStructs{:,3});
subStructs = {subStructs{:}}';

%Create data files
tic
for i = 1:length(noiseLevels)
    if i>length(noiseLevels) %sample the prior for one subStruct
        a = loadFieldData('subStructChoice',subStructs{i});
        % b = inversion(a,'priorOn',true);
    else
        a = loadFieldData('subStructChoice',subStructs{i});        
        b = inversion(a);
    end
end
toc

