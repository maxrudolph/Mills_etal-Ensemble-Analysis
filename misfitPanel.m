function outValues = misfitPanel(ewre2n,data,allModels,letter,noiseLevel)
allModels{1}.setWRE2N(data);
histogram(ewre2n,100,'EdgeAlpha',0);
outValues = zeros(1,size(allModels,2)+1);
hold on;
title(['\epsilon_n = ', noiseLevel])
yy = get(gca,'YLim');
[f,xi] = ksdensity(ewre2n);
[~,ind] = max(f);
normalizer = xi(ind);
for iPlot = 1:size(allModels,2)
    plot(allModels{iPlot}.wre2n*[1 1],yy,'LineStyle',allModels{iPlot}.lineStyle,...
        'Color',allModels{iPlot}.color,'LineWidth',1.0);
    outValues(iPlot) = allModels{iPlot}.wre2n/normalizer;
end
outValues(end) = normalizer;
set(gca,'FontSize',10);
xlabel('Weighted relative error');
text(0.95,0.9,letter,'units','normalized','FontSize',14);
bounds = get(gca,'XLim');
set(gca,'XLim',[bounds])
set(gca,'XScale','log')
end

