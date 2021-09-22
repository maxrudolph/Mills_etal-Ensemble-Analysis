clear;
close all;
filenames = {'3LayerA_0.2_02-Jul-2021.mat'};
titles = {'0.2'};
numEnsembles = length(filenames);

t = tiledlayout(3,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
figure1 = gcf();
figure1.Position(3:4) = [660 400];
set(gcf,'color','white');
load(['New folder/Ensemble_' filenames{1}],'results')

% Assumes the following order of plots:
%
C = [0 0 0;
    colororder()
    ];
line_widths = {1.5,1.5,1.5,1.5,1.5,1.5};
line_styles = {'-','-','-','-','-','-'};
ind = [1,2,3,4,5,6];
h=[];
i = 1;
load(['New folder/Analysis_' filenames{i}]);
load(['New folder/Ensemble_' filenames{i}],'results','data','forwardModel');
nexttile(5);
for j=1:length(allClusterSets{2}) % re-assign colors based on indexing into color order above
    allClusterSets{2}{j}.color = C(ind(j),:);
    allClusterSets{2}{j}.lineStyle = line_styles{j};
    %     if j>1
    %          allModels{j}.displayName = displayNames{j-1};
    %    end
end
importantNumbers = misfitPanel(ewre2n, results,data,forwardModel,allClusterSets{2},...
    2,titles{i},line_widths)
nexttile(1, [2 1])
modelSpacePanel(binCenters,numElements,allClusterSets{2},1,line_widths);
%     colormap(crameri('lajolla'));
colormap(flipud(gray));
legend('location','southeast')
lgd = legend('location','southeast');
ylabel('Depth (m)')
lgd.FontSize = 7;
c=colorbar();
c.Label.String = 'Probability (normalized)';

nexttile(6)
histogram(allPartitions{2}.indices);
xticks([1:5])
xlabel('Cluster #')
ylabel('Probability')
text(0.9,0.9,'D','units','normalized','FontSize',14);
set(gca,'YTick',[]);


nexttile(2, [2 1])

dataSpacePlot(xVals,yVals,allClusterSets{2},data);

text(0.9,0.9,'C','units','normalized','FontSize',14);


%{
%% Save the figure
set(gcf,'Visible','off');
set(gcf,'Renderer','painters');
exportgraphics(t,'test.eps');
set(gcf,'Renderer','opengl');
%}
set(gcf,'Visible','on');