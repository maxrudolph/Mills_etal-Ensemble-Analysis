clear;
close all;

% ensemble_files = {
% 'Ensemble_3LayerA_0.02.mat',
% 'Ensemble_3LayerA_0.1.mat',
% 'Ensemble_3LayerA_0.05.mat',
% 'Ensemble_3LayerA_0.01.mat',
% 'Ensemble_3LayerA_0.2.mat',
% 'Ensemble_3LayerA_0.mat'
% };
ensemble_files = {
  'Ensemble_Constable1984_1.mat'
% 'Ensemble_Constable1984_1_PRIOR_15-Nov-2023.mat'
};

delete(gcp('nocreate'));
parpool(3);
for i=1:length(ensemble_files)
  ensembleAnalysisMaster( ensemble_files{i},false);
end
