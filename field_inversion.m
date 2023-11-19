clear;
close all;

%function saveEnsemblesLoop
noiseLevels = [0,0,0,0]; % no additional noise.
subStructs = {'Constable1984_Wauchope','Constable1984_Renner','Constable1984_Wauchope','Constable1984_Renner'};
hierarchical = [true,true,false,false];%false,true,false];
prior = int32([2 2 2 2]);%1=flat 2=gaussian on rho
% nstruct = length(subStructs);
% nnoise = length(noiseLevels);
% noiseLevels = repmat(noiseLevels,1,nstruct)';
% subStructs = repmat(subStructs,nnoise,1);
% subStructs = vertcat(subStructs{:,1},subStructs{:,2},subStructs{:,3});
% subStructs = {subStructs{:}}';

% Create data files
delete(gcp('nocreate'));
parpool(4);
tic
parfor i = 1:length(noiseLevels)
    if i>length(noiseLevels) %sample the prior for one subStruct
        filename = loadFieldData('subStructChoice',subStructs{1});
        b = inversion(filename,'priorOn',true,'rhoPrior',prior(1));
    else
        filename = loadFieldData('subStructChoice',subStructs{i});       
        b = inversion(filename,'hierarchical',hierarchical(i),'rhoPrior',prior(i));
    end
end
toc

