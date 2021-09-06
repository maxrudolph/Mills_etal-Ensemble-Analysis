3filenames = {
         '3LayerA_0_02-Jul-2021.mat';
          '3LayerA_0.01_02-Jul-2021.mat';
    '3LayerA_0.02_02-Jul-2021.mat';
    '3LayerA_0.05_02-Jul-2021.mat';
    '3LayerA_0.1_02-Jul-2021.mat';
    '3LayerA_0.2_02-Jul-2021.mat'};
titles = {'0','0.01','0.02','0.05','0.1','0.2'}; 
tileLetter = 'A';
numEnsembles = length(filenames);
figure(1)
figure1 = tiledlayout(2,numEnsembles);
figure1.TileSpacing = 'compact';
figure1.Padding = 'compact';
set(gcf,'color','white');
figure(2)
figure2 = tiledlayout(2,numEnsembles/2);
figure2.TileSpacing = 'compact';
figure2.Padding = 'compact';
set(gcf,'color','white');

for ii = 1:numEnsembles
    disp('round')

load(['Ensemble_' filenames{ii}],'results')
    figure(1)
nexttile(ii);
histogram(results.ensembleNumLayers,'EdgeAlpha',0);
title(['\epsilon_n = ',titles(ii)]),xlabel('# of layers'),ylabel('# of solutions')
%set(gca,'YScale','log')
nexttile(ii+numEnsembles)
histogram(results.ensembleNumLayers,'EdgeAlpha',0);
title(['\epsilon_n = ',titles(ii)]),xlabel('# of layers'),ylabel('# of solutions')
set(gca,'YScale','log');

figure(2)
nexttile(ii)
histogram(log10(results.ensembleVars),'EdgeAlpha',0);
title('Ensemble Variance'),xlabel('log(\sigma^2)'),ylabel('# of solutions')
end
