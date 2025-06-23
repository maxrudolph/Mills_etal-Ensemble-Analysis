clear;
close all;

ensemble_files = {    
'Ensemble_3LayerA__hierarchical-1_rhoPrior-1_0.02.mat',
'Ensemble_3LayerA__hierarchical-1_rhoPrior-1_0.05.mat',
'Ensemble_3LayerA__hierarchical-1_rhoPrior-1_0.1.mat',
'Ensemble_3LayerA__hierarchical-1_rhoPrior-1_0.mat'
};
exact_known = true;

% delete(gcp('nocreate'));
% parpool();

for i=1:length(ensemble_files)
    ensembleAnalysisMaster_noClustering( ensemble_files{i},exact_known);
end
