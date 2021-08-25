clear
filenames = {
%     '3LayerA_0_02-Jul-2021.mat';
%      '3LayerA_0.01_02-Jul-2021.mat';
       '3LayerA_0.02_02-Jul-2021.mat';
    '3LayerA_0.05_02-Jul-2021.mat';
    '3LayerA_0.1_02-Jul-2021.mat'};
%    '3LayerA_0.2_02-Jul-2021.mat'};
titles = {'0.02','0.05','0.1'}; %'0','0.01',
theAlphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
numEnsembles = length(filenames);

figure1 = tiledlayout(3,numEnsembles);
figure1.TileSpacing = 'compact';
figure1.Padding = 'compact';
set(gcf,'color','white');
load(['Ensemble_' filenames{1}],'results')

for i = 1:numEnsembles
    disp('round')
    load(['Analysis_' filenames{i}]);
    load(['Ensemble_' filenames{i}],'results','data','forwardModel');
    nexttile(i)
    importantNumbers = misfitPanel(ewre2n, results,data,forwardModel,allModels,...
        theAlphabet(2*i-1),titles{i})
    nexttile(i+numEnsembles,[2 1])
    modelSpacePanel(binCenters,numElements,allModels,theAlphabet(2*i));
    if i == 1
        legend('location','southeast')
lgd = legend('location','southeast');
lgd.FontSize = 7;
    end
end


nexttile(1)
ylabel('Solutions');
nexttile(numEnsembles+1);
ylabel('Depth (m)');
nexttile(numEnsembles*2)
c=colorbar();
c.Label.String = 'Relative normalized probability';