types = {'1LayerA','3LayerA','4LayerA'};
noiseLevels = {'0','0.01','0.02','0.05','0.1','0.2','prior'};
addpath(genpath(fileparts(mfilename('fullpath'))))

for i = 1:4
    mkdir(['Allfigure',num2str(i)]);
end

for i = 1:length(types)
    for j = 1:length(noiseLevels)

        filename = ['/work/cdmills/Senior-Thesis-Stuff/Ensembles_07022021/' 'Ensemble_',types{i},'_',noiseLevels{j},'_02-Jul-2021.mat']

        ensembleAnalysisA(filename,true)
        

    end
end
disp('All done');