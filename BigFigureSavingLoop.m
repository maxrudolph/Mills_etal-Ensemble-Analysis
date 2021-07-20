types = {'1LayerA','3LayerA','4LayerA'};
noiseLevels = {'0','0.01','0.02','0.05','0.1','0.2'};
addpath(genpath(fileparts(mfilename('fullpath'))))

parfor i = 1:length(types)
    for j = 1:length(noiseLevels)

        filename = ['/work/cdmills/Senior-Thesis-stuff/Ensembles_07022021/' 'Ensemble_',types{i},'_',noiseLevels{j},'_02-Jul-2021.mat']

        ensembleAnalysisMaster(filename)
        
    end
end

filename = ['/work/cdmills/Senior-Thesis-stuff/Ensembles_07022021/'...
    'Ensemble_1LayerA_0_PRIOR_02-Jul-2021.mat'];
ensembleAnalysisMaster(filename);

disp('All done');