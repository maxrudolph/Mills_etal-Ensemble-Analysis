clear;
close all;

%function saveEnsemblesLoop
noiseLevels = [0]; % no additional noise.
% subStructs = {'Constable1984_Renner','Constable1984_Renner','Constable1984_Renner','Constable1984_Renner'};
subStructs = {'Constable1984_Wauchope','Constable1984_Wauchope','Constable1984_Wauchope','Constable1984_Wauchope'};
nstruct = length(subStructs);
nnoise = length(noiseLevels);
% noiseLevels = repmat(noiseLevels,1,nstruct)';
% subStructs = repmat(subStructs,nnoise,1);
% subStructs = vertcat(subStructs{:,1},subStructs{:,2},subStructs{:,3});
% subStructs = {subStructs{:}}';
priors = int64([1,2,1,2]);
sample_prior = [true, true, false, false];

delete(gcp('nocreate'));
parpool(4);
%Create data files
tic
parfor i = 1:4%1:length(noiseLevels)
   
        a = loadFieldData('subStructChoice',subStructs{i});        
        b = inversion(a,'priorOn',sample_prior(i),'priorChoice',priors(i));
   
end
toc
