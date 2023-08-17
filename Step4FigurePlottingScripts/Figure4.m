clear;
close all;

% file_prefix = '../Ensembles_09132021/';
% file_prefix = '../Ensembles_02082023/';
file_prefix = './';

% file_prefix = '~/Box/Davis/Students/Chris Mills/MCMC Box Shared Folder/Ensembles/Ensembles_09132021/';
filenames = {
    %'3LayerA_0_02-Jul-2021.mat';
    %'3LayerA_0.01_02-Jul-2021.mat';
    '3LayerA_0.02.mat';
    '3LayerA_0.05.mat';
    '3LayerA_0.1.mat';
    %     '3LayerA_0.2_02-Jul-2021.mat'
    };
titles = {'0.02','0.05','0.1','0.2','-'};
numEnsembles = length(filenames);

t = tiledlayout(3,numEnsembles);
t.TileSpacing = 'compact';
t.Padding = 'compact';
figure1 = gcf();
figure1.Position(3:4) = [660 400];
set(gcf,'color','white');
% load([file_prefix 'Ensemble_' filenames{1}],'results')

% Assumes the following order of plots:
% 'Exact solution', 'MS Mean','MS Median','MS Max Likelihood'.'DS Best Fit','DS Median'
C = [0 0 0;
    colororder()
    ];

line_widths = {1.5,1.5,1.5,1.5,1.5,1.5};
line_styles = {'-','-','--','-','--','-'};
ind = [1,2,2,4,4,2,3,4];
displayNames = {'k-Medoids Euclidian 1','k-Medoids Euclidean 2','k-Medoids Manhattan 1','k-Medoids Manhattan 2',' ',' ',' '};

h=[];
for i = 1:numEnsembles
    load([file_prefix 'Analysis_' filenames{i}]);
    load([file_prefix 'Ensemble_' filenames{i}],...
        'results','data','forwardModel');
    for j=1:length(allClusterSets)
        clusterset_weighted_errors = cellfun( @(x) x.wre2n,allClusterSets{j} );
        [~,ind1] = sort(clusterset_weighted_errors(2:end)); %sort weighted errors in 2nd position through end - 1st position is true solution, 2nd through end are clustering results.
        allClusterSets{j}(2:end) = allClusterSets{j}(ind1+1);
    end
    allModels = {allClusterSets{3}{:},allClusterSets{4}{2:end}};
    nexttile(i)
    for j=1:length(allModels) % re-assign colors based on indexing into color order above
        allModels{j}.color = C(ind(j),:);
        allModels{j}.lineStyle = line_styles{j};
        if j>1
            allModels{j}.displayName = displayNames{j-1};
        end
    end
    titles{i}
    importantNumbers = misfitPanel(ewre2n, results,data,forwardModel,allModels,...
        2*i-1,titles{i},line_widths);
    cellfun( @(x) x.displayName,allModels,'UniformOutput',false)
    round(importantNumbers,2)
    nexttile(i+numEnsembles,[2 1])
    modelSpacePanel(binCenters,numElements,allModels,2*i,line_widths);
    %     colormap(crameri('lajolla'));
    colormap(flipud(gray));
    if i == 2
        legend('location','southeast')
        lgd = legend('location','southeast');
        lgd.FontSize = 7;
    end
end
nexttile(1); xticks([.01 .02 .03 .04 .05 .1]);
nexttile(2); xticks([.05 0.07 .1 .2 .3]);
nexttile(3); xticks(0.1:0.1:0.4);

nexttile(1)
ylabel('Probability');
nexttile(numEnsembles+1);
ylabel('Depth (m)');
nexttile(numEnsembles*2)
c=colorbar();
c.Label.String = 'Probability (normalized)';


%% Save the figure
figure(figure1);
set(gcf,'Visible','off');
set(gcf,'Renderer','painters');
exportgraphics(t,'Figure4.eps');
set(gcf,'Renderer','opengl');

set(gcf,'Visible','on');