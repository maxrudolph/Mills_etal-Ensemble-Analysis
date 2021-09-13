types = {'3LayerA'};
noiseLevels = {'0.1','0.05','0.02','0.2','0.01','0'};
addpath(genpath(fileparts(mfilename('fullpath'))))


for i = 1:length(types)
    parfor j = 1:length(noiseLevels)

        filename = [...
            '/work/cdmills/Senior-Thesis-stuff/Ensembles_07022021/Ensemble_',...
            types{i},'_',noiseLevels{j},'_02-Jul-2021.mat']

        ensembleAnalysisMaster(filename);
        
    end
end

%filename = ['_1LayerA_0_PRIOR_02-Jul-2021.mat'];
%figurePlotting(filename,true);

disp('All done');