types = {'1LayerA','3LayerA','4LayerA'};
noiseLevels = {'0.01','0.02','0.05','0.1','0.2'};
addpath(genpath(fileparts(mfilename('fullpath'))))

for i = 1:6
    mkdir(['Allfigure' num2str(i)]);
end

parfor i = 1:length(types)
    for j = 1:length(noiseLevels)

        filename = ['_',types{i},'_',noiseLevels{j},'_02-Jul-2021']

        figurePlotting(filename,true)
        
    end
end

filename = ['_1LayerA_0_PRIOR_02-Jul-2021.mat'];
figurePlotting(true);

disp('All done');