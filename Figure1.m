clear
filenames = {'3LayerA_0.02_02-Jul-2021.mat';
    '3LayerA_0.05_02-Jul-2021.mat';
    '3LayerA_0.1_02-Jul-2021.mat';
    '3LayerA_0.2_02-Jul-2021.mat'};
titles = {'0.02','0.05','0.1','0.2'};
theAlphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
numEnsembles = length(filenames);

figure1 = tiledlayout(3,numEnsembles);
figure1.TileSpacing = 'compact';
figure1.Padding = 'compact';
set(gcf,'color','white');

for i = 1:numEnsembles
    disp('round')
    load(['Analysis_' filenames{i}]);
    load(['Ensemble_' filenames{i}],'results','data','forwardModel');
    nexttile(i)
    misfitPanel(results,data,forwardModel,allClusterSets{1},theAlphabet(2*i-1),...
        titles{i});
    nexttile(i+numEnsembles,[2 1])
    modelSpacePanel(binCenters,numElements,allClusterSets{1},theAlphabet(2*i));
end


nexttile(1)
ylabel('Solutions');
nexttile(numEnsembles+1);
ylabel('Depth (m)');
nexttile(numEnsembles*2)
c=colorbar();
c.Label.String = 'Relative normalized probability';