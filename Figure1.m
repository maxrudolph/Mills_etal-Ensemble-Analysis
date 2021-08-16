figure1 = tiledlayout(3,3)

filenames = {'Analysis_3LayerA_0.11_15-Jul-2021.mat';
    'Analysis_1LayerA_0.05_02-Jul-2021.mat';
    'Analysis_1LayerA_0_02-Jul-2021.mat'};
load('Ensemble_3LayerA_0.11_15-Jul-2021.mat')
xx = zeros(3,2);

for i = 1:3
    disp('round')
    load(filenames{i});
    nexttile
    histogram(results.ensembleMisfits,100,'EdgeAlpha',0);
    hold on;
    yy=get(gca,'YLim');
    xx(i,:) = get(gca,'XLim')
    for iPlot = 2:size(allModels,2)
    plot(allModels{iPlot}.misfit*[1 1],yy,'LineStyle',...
        allModels{iPlot}.lineStyle,'Color',allModels{iPlot}.color,...
        'LineWidth',1.0);
    end
    set(gca,'FontSize',10);
    if i ==2
        xlabel('Misfit (\Omega-m)');
    end
    if i ==1
        ylabel('# of solutions')
    end
    
    nexttile(i+3,[2 1]);
    
    xdata = 10.^binCenters{1};
    mask = xdata >= min(allModels{1}.rhos)/100 & xdata <= max(allModels{1}.rhos)*100;
    % xlim([min(allModels{1}.rhos)/100,max(allModels{1}.rhos)*100]);
    xticks = 10.^(floor(log10(min(xdata))):2:ceil(log10(max(xdata))));
    p = pcolor(xdata(mask),10.^binCenters{2},numElements(mask,:)'); shading flat;

    colormap(flipud(bone))
    set(gca,'XScale','log','YScale','log','ColorScale','log');
    hold on
for iPlot2 = 1:size(allModels,2)
    plot(allModels{iPlot2}.rhos,allModels{iPlot2}.depths,'LineStyle',...
        allModels{iPlot2}.lineStyle,'Color',allModels{iPlot2}.color,'DisplayName',...
        allModels{iPlot2}.displayName,'LineWidth',1.0);
end
ylim = get(gca,'YLim');
set(gca,'YLim',[ylim(1) ylim(2)]);
if i==3
    c=colorbar();
    c.Label.String = 'Number of solutions';
end
set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%legend();
set(gca,'YDir','reverse','FontSize',10,'Box','on','Layer','Top');
k = find(sum(numElements')>0);

ax = gca;
ax.XTickMode = 'manual';
ax.XTick = xticks;
if i==2
    xlabel('Resistivity (\Omega-m)');
end
if i ==1
    ylabel('Depth (m)');
end
%lgd = legend('location','southeast');
%lgd.FontSize = 7;
%text(0.9,0.95,'C','units','normalized','FontSize',14)
    
end

for i = 1:3
    nexttile(i)
    set(gca,'XLim',[min(xx(:,1)) max(xx(:,2)) ]);
end
