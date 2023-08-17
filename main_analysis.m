clear;
close all;

ensemble_files = {
'Ensemble_3LayerA_0.02.mat',
'Ensemble_3LayerA_0.1.mat',
'Ensemble_3LayerA_0.05.mat',
'Ensemble_3LayerA_0.01.mat',
'Ensemble_3LayerA_0.2.mat',
'Ensemble_3LayerA_0.mat'
};

delete(gcp('nocreate'));
for i=1:length(ensemble_files)
  ensembleAnalysisMaster( ensemble_files{i} );  
end
