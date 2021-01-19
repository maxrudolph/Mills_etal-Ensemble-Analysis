types = {'1LayerA','3LayerA','4LayerA'};
noiseLevels = {'0','0.01','0.02','0.05','0.1','0.2','prior'};

for i = 1:8
    mkdir(['Allfigure',num2str(i)]);
end

for i = 1:length(types)
    for j = 1:length(noiseLevels)
        filename = ['Ensemble_',types{i},'_',noiseLevels{j},'_02-Dec-2020.mat']
        ensembleAnalysis3
    end
end
disp('All done');