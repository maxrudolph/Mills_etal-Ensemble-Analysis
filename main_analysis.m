clear;
close all;

ensemble_files = {
    'Ensemble_Constable1984_Renner__hierarchical-1_rhoPrior-1_1_22-Feb-2024.mat',
    'Ensemble_Constable1984_Renner__hierarchical-1_rhoPrior-2_1_22-Feb-2024.mat'  
};

delete(gcp('nocreate'));
parpool(3);

for i=1:length(ensemble_files)
    ensembleAnalysisMaster( ensemble_files{i},true);
end
