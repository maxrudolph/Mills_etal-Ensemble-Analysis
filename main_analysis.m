clear;
close all;

ensemble_files = {    
%'Ensemble_3LayerA__hierarchical-1_rhoPrior-1_0.02.mat',
%'Ensemble_3LayerA__hierarchical-1_rhoPrior-1_0.05.mat',
% 'Ensemble_3LayerA__hierarchical-1_rhoPrior-1_0.1.mat',
%'Ensemble_3LayerA__hierarchical-1_rhoPrior-1_0.mat'
'Ensemble_Constable1984_Wauchope__hierarchical-1_rhoPrior-2_1_PRIOR',
'Ensemble_Constable1984_Wauchope__hierarchical-1_rhoPrior-1_1_PRIOR'
};
exact_known = false;

% delete(gcp('nocreate'));
% parpool();

for i=1:length(ensemble_files)
    ensembleAnalysisMaster_noClustering( ensemble_files{i},exact_known);
end
