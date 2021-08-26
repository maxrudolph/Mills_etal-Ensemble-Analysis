clear;
close all;
filenames = {
    %     '3LayerA_0_02-Jul-2021.mat';
    %      '3LayerA_0.01_02-Jul-2021.mat';
    '3LayerA_0.02_02-Jul-2021.mat';
    '3LayerA_0.05_02-Jul-2021.mat';
    '3LayerA_0.1_02-Jul-2021.mat'};
%    '3LayerA_0.2_02-Jul-2021.mat'};
titles = {'0.02','0.05','0.1'}; %'0','0.01',
numEnsembles = length(filenames);

t = tiledlayout(3,numEnsembles);
t.TileSpacing = 'compact';
t.Padding = 'compact';
figure1 = gcf();
figure1.Position(3:4) = [660 400];
set(gcf,'color','white');
load(['Ensemble_' filenames{1}],'results')

% Assumes the following order of plots:
% 'Exact solution', 'MS Mean','MS Median','MS Max Likelihood'.'DS Best Fit','DS Median'
C = [0 0 0;
   colororder()     
    ]
line_widths = {1.5,1.5,1.5,2,1.5,1.5};
line_styles = {'-','--','--','--','-','-'};
ind = [1,3,4,7,2,4];


h=[];
for i = 1:numEnsembles    
    load(['Analysis_' filenames{i}]);
    load(['Ensemble_' filenames{i}],'results','data','forwardModel');
    nexttile(i)
    for j=1:length(allModels) % re-assign colors based on indexing into color order above
        allModels{j}.color = C(ind(j),:);
        allModels{j}.lineStyle = line_styles{j};
    end
    importantNumbers = misfitPanel(ewre2n, results,data,forwardModel,allModels,...
        2*i-1,titles{i},line_widths)
    nexttile(i+numEnsembles,[2 1])
    modelSpacePanel(binCenters,numElements,allModels,2*i,line_widths);
%     colormap(crameri('lajolla'));
    colormap(flipud(gray));
    if i == 1
        legend('location','southeast')
        lgd = legend('location','southeast');
        lgd.FontSize = 7;
    end
end
nexttile(1); xticks([.01 .02 .03 .04]);
nexttile(2); xticks([.05 .1 .2]);
nexttile(3); xticks(0.1:0.1:0.4);

nexttile(1)
ylabel('Probability');
nexttile(numEnsembles+1);
ylabel('Depth (m)');
nexttile(numEnsembles*2)
c=colorbar();
c.Label.String = 'Probability (normalized)';

%% Save the figure
set(gcf,'Visible','off');
set(gcf,'Renderer','painters');
exportgraphics(t,'test.eps');
set(gcf,'Renderer','opengl');
set(gcf,'Visible','on');