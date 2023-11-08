function outValues = misfitPanel(ewre2n,results,data,forwardModel,allModels,panel_number,noiseLevel,line_widths)
% note - outValues is a list of the weighted residual in the 2-norm for
% each model in allModels. The last entry in outValues is the mode of the
% weighted relative error (used for normalization)

histogram(ewre2n,100,'EdgeAlpha',0,'FaceColor',0.65*[1 1 1]);
if ~isempty(allModels) && ~isempty(allModels{1})
    allModels{1}.setWRE2N(data);
end

outValues = zeros(1,size(allModels,2)+1);
hold on;
title(['\epsilon_n = ', noiseLevel])
yy = get(gca,'YLim');
[f,xi] = ksdensity(ewre2n);
[~,ind] = max(f);
normalizer = xi(ind);
for iPlot = 1:size(allModels,2)
    if isfield(allModels{iPlot},'wre2n')
        plot(allModels{iPlot}.wre2n*[1 1],yy,'LineStyle',allModels{iPlot}.lineStyle,...
            'Color',allModels{iPlot}.color,'LineWidth',line_widths{iPlot});
        outValues(iPlot) = allModels{iPlot}.wre2n/normalizer;
    end
end
outValues(end) = normalizer;
set(gca,'FontSize',10);
xlabel('Weighted relative error');
text(0.9,0.9,char(64+panel_number),'units','normalized','FontSize',14);
% bounds = get(gca,'XLim');
% set(gca,'XLim',[bounds])
set(gca,'XScale','log');
bounds = get(gca,'XLim');
set(gca,'XLim',[min(bounds*0.95) max(bounds*1.05)]);
% ticks = get(gca,'XTick');

set(gca,'YTick',[]);
end

