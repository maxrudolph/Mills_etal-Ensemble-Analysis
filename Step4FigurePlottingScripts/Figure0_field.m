clear;
close all;
addpath Step3AnalysisScripts
addpath Step2InversionScripts
addpath Step4FigurePlottingScripts

% file_prefix = '~/Box/Davis/Students/Chris Mills/MCMC Box Shared Folder/Ensembles/Ensembles_09102021/';
file_prefix = './'
% file_prefix='../Ensembles_02082023/'
% file_prefix = '../Ensembles_09132021/'
filenames = {
    'Constable1984_1_PRIOR_09-Nov-2023.mat'
    %'3LayerA_0_02-Jul-2021.mat';
%     '3LayerA_0.05.mat';
     % '3LayerA_0.02.mat';
     % '3LayerA_0.05.mat';
     % '3LayerA_0.1.mat';
%     '3LayerA_0.2_02-Jul-2021.mat'
};
titles = {'0.02','0.05','0.1'};
% titles={'0.05'};
numEnsembles = length(filenames);

t = tiledlayout(5,numEnsembles);
t.TileSpacing = 'compact';
t.Padding = 'compact';
figure1 = gcf();
figure1.Position(3:4) = [660 400];
set(gcf,'color','white');
load([file_prefix 'Ensemble_' filenames{1}],'results')

% Assumes the following order of plots:
% 'Exact solution', 'MS Mean','MS Median','MS Max Likelihood'.'DS Best Fit','DS Median'
C = [0 0 0;
    colororder()
    ]
line_widths = {1.5,1.5,1.5,1.5,1.5,1.5};
line_styles = {'-','--','--','--','-','-'};
ind = [1,3,4,7,2,4];

h=[];
for i = 1:numEnsembles
    figure(figure1);
    % load([file_prefix 'Analysis_' filenames{i}]);
    load([file_prefix 'Ensemble_' filenames{i}],'results','data','forwardModel');
    nexttile(i)
    for j=1:length(allModels) % re-assign colors based on indexing into color order above
        allModels{j}.color = C(ind(j),:);
        allModels{j}.lineStyle = line_styles{j};
        %     if j>1
        %          allModels{j}.displayName = displayNames{j-1};
        %    end
    end
    titles{i}
    importantNumbers = misfitPanel(ewre2n, results,data,forwardModel,[],...
        2*i-1,titles{i},line_widths)
    % histogram of number of layers
    nexttile(i+numEnsembles)
    histogram(results.ensembleNumLayers,'BinEdges',0.5:10.5,'FaceColor',0.65*[1 1 1]);
    set(gca,'YTick',[]);
    text(0.90,0.95,char(64+2*i),'units','normalized','FontSize',14);
    set(gca,'XTick',1:10);
    
    % histograms of variance
    nexttile(i+2*numEnsembles);
    histogram(results.ensembleVars,'FaceColor',0.65*[1 1 1],'EdgeColor','none');
    hold on
    plot(str2num(titles{i})^2*[1 1],get(gca,'YLim'),'r')
    xlabel('Noise Hyperparameter')
    
    nexttile(i+3*numEnsembles,[2 1])
    % model space pdf
    modelSpacePanel(binCenters,numElements,[],3*i,line_widths);
    set(gca,'ColorScale','linear');
    %     colormap(crameri('lajolla'));
    colormap(flipud(gray));

    % diagnostic plot for misfit
    figure();
    n3mask = results.ensembleNumLayers==3;
    n4mask = results.ensembleNumLayers==4;
    n5mask = results.ensembleNumLayers==5;
    histogram(results.ensembleMisfits(n3mask));
    hold on;
    histogram(results.ensembleMisfits(n4mask),'EdgeColor','none','FaceColor','r','FaceAlpha',0.5);
    histogram(results.ensembleMisfits(n5mask),'EdgeColor','none','FaceColor','g','FaceAlpha',0.5);
    
    set(gca,'XScale','log')
    title(titles{i});
    
    
end
% nexttile(1); xticks([.01 .02 .03 .04]);
% nexttile(2); xticks([.05 .1 .2]);
% nexttile(3); xticks(0.1:0.1:0.4);

nexttile(1)
ylabel('Probability');
nexttile(numEnsembles+1);
ylabel('Frequency (-)');
nexttile(3*numEnsembles+1);
ylabel('Depth (m)');
nexttile(numEnsembles*4)
c=colorbar();
c.Label.String = 'Probability (normalized)';

%{
%% Save the figure
set(gcf,'Visible','off');
set(gcf,'Renderer','painters');
exportgraphics(t,'test.eps');
set(gcf,'Renderer','opengl');
%}
set(gcf,'Visible','on');