function misfitPanel(results,data,forwardModel,allModels,letter,noiseLevel)
ensembleChiSquareds = zeros(1,length(results.ensembleMisfits));
Cdi = inv(data.Cd);
denominator = data.y'*Cdi*data.y;
for iPlot = 1:length(results.ensembleMisfits)
    residual = forwardModel(results.ensembleDepths(:,iPlot),...
        results.ensembleRhos(:,iPlot),data.lambda) - data.y;
    ensembleChiSquareds(iPlot) = 1 - ((residual'*Cdi*residual)/...
        (denominator));
end
histogram(ensembleChiSquareds,100,'EdgeAlpha',0);
hold on;
title(noiseLevel)
yy = get(gca,'YLim');
for iPlot = 1:size(allModels,2)
    residual = allModels{iPlot}.y - data.y;
    chiSquared = 1 - ((residual'*Cdi*residual)/denominator);
    plot(chiSquared*[1 1],yy,'LineStyle',allModels{iPlot}.lineStyle,...
        'Color',allModels{iPlot}.color,'LineWidth',1.0);
end
set(gca,'FontSize',10);
xlabel('\chi^2 distribution');
text(0.05,0.9,letter,'units','normalized','FontSize',14);
bounds = get(gca,'XLim');
set(gca,'XLim',[bounds(1),1])
end

