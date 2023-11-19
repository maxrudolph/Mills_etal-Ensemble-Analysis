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
    'Ensemble_Constable1984_Wauchope__hierarchical-1_rhoPrior-2_1_18-Nov-2023.mat',
    'Ensemble_Constable1984_Wauchope__hierarchical-0_rhoPrior-2_1_18-Nov-2023.mat',
    'Ensemble_Constable1984_Renner__hierarchical-1_rhoPrior-2_1_18-Nov-2023.mat',
    'Ensemble_Constable1984_Renner__hierarchical-0_rhoPrior-2_1_18-Nov-2023.mat'
};

delete(gcp('nocreate'));
parpool(3);
for i=1:length(ensemble_files)
  ensembleAnalysisMaster_noClustering( ensemble_files{i},false);
end
