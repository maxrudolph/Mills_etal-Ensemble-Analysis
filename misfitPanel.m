function ewre2n = misfitPanel(ewre2n,results,data,forwardModel,allModels,letter,noiseLevel)
histogram(ewre2n,100,'EdgeAlpha',0);
hold on;
title(['\epsilon_n = ', noiseLevel])
yy = get(gca,'YLim');
normalizer = mean(ewre2n);
for iPlot = 1:size(allModels,2)
    plot(allModels{iPlot}.wre2n*[1 1],yy,'LineStyle',allModels{iPlot}.lineStyle,...
        'Color',allModels{iPlot}.color,'LineWidth',1.0);
end
set(gca,'FontSize',10);
xlabel('Weighted relative error');
text(0.95,0.9,letter,'units','normalized','FontSize',14);
bounds = get(gca,'XLim');
set(gca,'XLim',[bounds])
set(gca,'XScale','log')
end

