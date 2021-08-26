function outValues = misfitPanel(ewre2n,results,data,forwardModel,allModels,panel_number,noiseLevel,line_widths)
histogram(ewre2n,100,'EdgeAlpha',0,'FaceColor',0.65*[1 1 1]);

outValues = zeros(1,size(allModels,2)+1);
hold on;
title(['\epsilon_n = ', noiseLevel])
yy = get(gca,'YLim');
[f,xi] = ksdensity(ewre2n);
[~,ind] = max(f);
normalizer = xi(ind);
for iPlot = 1:size(allModels,2)
    plot(allModels{iPlot}.wre2n*[1 1],yy,'LineStyle',allModels{iPlot}.lineStyle,...
        'Color',allModels{iPlot}.color,'LineWidth',line_widths{iPlot});
    outValues(iPlot) = allModels{iPlot}.wre2n/normalizer;
end
outValues(end) = normalizer;
set(gca,'FontSize',10);
xlabel('Weighted relative error');
text(0.9,0.9,char(64+panel_number),'units','normalized','FontSize',14);
bounds = get(gca,'XLim');
set(gca,'XLim',[bounds])
set(gca,'XScale','log')
end

