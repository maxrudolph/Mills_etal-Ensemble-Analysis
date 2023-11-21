clear;
close all;
addpath Step1DataGenerationScripts
addpath Step3AnalysisScripts
addpath Step2InversionScripts
addpath Step4FigurePlottingScripts

file_prefix = './'


filenames = {      
'Constable1984_Wauchope__hierarchical-0_rhoPrior-2_1.mat',
'Constable1984_Wauchope__hierarchical-1_rhoPrior-2_1'

    };
titles = {'wauchope-1','wauchope-0','renner-1','renner-0'};
% titles={'0.05'};
numEnsembles = length(filenames);

t = tiledlayout(6,numEnsembles);
t.TileSpacing = 'compact';
t.Padding = 'compact';
figure1 = gcf();
figure1.Position(3:4) = [275*2 720]*get(groot,'ScreenPixelsPerInch')/72;
set(gcf,'color','white');
% load([file_prefix 'Ensemble_' filenames{1}],'results')

exact_color = 'r';
observations_color = [80 157 191]/255;

% Assumes the following order of plots:
% 'Exact solution', 'MS Mean','MS Median','MS Max Likelihood'.'DS Best Fit','DS Median'
C = [0 0 0;
    colororder()
    ]
line_widths = {0.5,1.5,1.5,1.5,1.5,1.5};
line_styles = {'-','--','--','--','-','-'};
ind = [1,3,4,7,2,4];

h=[];
for i = 1:numEnsembles    
    load([file_prefix 'Analysis_' filenames{i}]);
    load([file_prefix 'Ensemble_' filenames{i}],'results','data','forwardModel');
    for j=1:length(allModels) % re-assign colors based on indexing into color order above
        allModels{j}.color = C(ind(j),:);
        allModels{j}.lineStyle = line_styles{j};
    end
    disp(['now on ' titles{i}]);
    %% Misfit - second row
    figure(figure1);
    nexttile(numEnsembles+i)
    histogram(ewre2n,250,'FaceColor',0.65*[1 1 1],'EdgeColor','none');
    text(0.90,0.90,char(64+5*(i-1)+2),'units','normalized','FontSize',14);
    set(gca,'YTick',[]);
    set(gca,'XLim',[0.0 0.2])
    hold on;
    % plot(allModels{1}.wre2n*[1 1],get(gca,'YLim'),'Color',exact_color); % add line for wre2n
    xlabel('Weighted Relative Error');
    %     importantNumbers = misfitPanel(ewre2n, results,data,forwardModel,[],...
    %         i,titles{i},line_widths)
    
    %% make 2D histogram of forward model predictions
    % generate high-resolution data space values
    nxplot = 1001;
    ntot = size(results.ensembleG,2);
    xvals = logspace(0,log10(2e5),nxplot);
    lamplot = makeLambda(xvals,11);
    Gplot = zeros(nxplot,ntot);

    for j=1:size(results.ensembleG,2)
        if ~mod(j,10000)
            disp([num2str(j/ntot*100) ' percent done'])
        end
        nlayer = results.ensembleNumLayers(j);
        Gplot(:,j) = forwardModel(results.ensembleDepths(1:nlayer,j),results.ensembleRhos(1:nlayer,j),lamplot);
    end
    
    %% bin data-space values and create a histogram (first row)
    figure(figure1);
    nexttile(i)
    
    yvals = logspace(0,4,1001);
    N = zeros(length(yvals)-1,nxplot);
    xplot = repmat(xvals',[1 ntot]);
    for j=1:nxplot
        N(:,j) = histcounts(Gplot(j,:),yvals);
    end
    %     xc = 10.^((log10(xbin(1:end-1)) + log10(xbin(2:end)))/2);
    yc = 10.^((log10(yvals(1:end-1)) + log10(yvals(2:end)))/2);
    imagesc(xvals,yc,N); hold on;
    set(gca,'YDir','normal');
    set(gca,'YLim',[5 max(data.y)*2]);
    colormap(flipud(gray));
    set(gca,'XScale','log','YScale','log');
    hold on
    plot(data.x,data.y,'.','MarkerFaceColor',observations_color,'MarkerEdgeColor',observations_color)
    set(gca,'XLim',[1 max(data.x)*2]);
    set(gca,'XTick',10.^[0 1 2 3 4 5]);
    set(gca,'YTick',10.^[0 1 2 3 4 5]);
    xlabel('Spacing (m)');
    text(0.90,0.90,char(64+(i-1)*5+1),'units','normalized','FontSize',14);
    title(['\epsilon_n = ' titles{i}]);

    %
    %% histogram of number of layers
    %
    nexttile(i+2*numEnsembles)
    histogram(results.ensembleNumLayers,'BinEdges',0.5:30.5,'FaceColor',0.65*[1 1 1]);
    set(gca,'YTick',[]);
    text(0.90,0.90,char(64+5*(i-1)+3),'units','normalized','FontSize',14);
    set(gca,'XTick',1:2:30);
    hold on
    %plot([3 3],get(gca,'Ylim'),'r');
    xlabel('Number of Layers');
    
    %
    %% histogram of noise hyperparameter
    %
    nexttile(i+3*numEnsembles);
    histogram(results.ensembleVars,'FaceColor',0.65*[1 1 1],'EdgeColor','none');
    hold on
    %plot(str2num(titles{i})^2*[1 1],get(gca,'YLim'),'Color',observations_color,'LineWidth',1)
    text(0.90,0.90,char(64+5*(i-1)+4),'units','normalized','FontSize',14);
    set(gca,'YTick',[]);
    % set(gca,'XLim',[0.0 0.2]);
    set(gca,'XLim',[1e-1 1e1],'XScale','log','XTick',[1e-4 1e-3 1e-2 1e-1 1e0 1e1]);
    xlabel('Noise Hyperparameter')
    
    %% model space pdf
    nexttile(i+4*numEnsembles,[2 1])
    allModels{1}.color = 'r';
    allModels{1}.lineWidth = 0.5;
    modelSpacePanel(binCenters,numElements,{allModels{1}},5*i,line_widths);
    set(gca,'ColorScale','linear');
    %     colormap(crameri('lajolla'));
    colormap(flipud(gray));
    set(gca,'XLim',[1e-1 1e7]);
    
    
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
    drawnow();
    
end
% nexttile(1); xticks([.01 .02 .03 .04]);
% nexttile(2); xticks([.05 .1 .2]);
% nexttile(3); xticks(0.1:0.1:0.4);
figure(figure1);
nexttile(1)
ylabel('A.R. (\Omega-m)');
nexttile(numEnsembles+1);
ylabel('Frequency');
nexttile(2*numEnsembles+1);
ylabel('Frequency');
nexttile(3*numEnsembles+1);
ylabel('Frequency');
nexttile(4*numEnsembles+1);
ylabel('Depth (m)');
nexttile(numEnsembles*5)
c=colorbar();
c.Label.String = 'Probability (normalized)';


%% Save the figure
figure(figure1);
disp('Saving...');
set(figure1,'Visible','off');
set(figure1,'Renderer','painters');
exportgraphics(t,'Figure1_field.eps');
set(figure1,'Renderer','opengl');
set(figure1,'Visible','on');
