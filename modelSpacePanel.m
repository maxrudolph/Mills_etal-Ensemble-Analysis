function h = modelSpacePanel(binCenters,numElements,allModels,panel_number,line_widths)

xdata = 10.^binCenters{1};
mask = xdata >= min(allModels{1}.rhos)/100 & xdata <= max(allModels{1}.rhos)*100;
% Pseudocolor plot of model-space probability density
p = pcolor(xdata(mask),10.^binCenters{2},((1/max(max(numElements)))*(numElements(mask,:)')));
shading flat;
cmap = colormap(flipud(bone));
set(gca,'XScale','log','YScale','log','ColorScale','log');
hold on

% add curves for model space solutions
for iPlot = 1:length(allModels) % iterate backwards to plot exact solution last.
    plot(allModels{iPlot}.rhos,allModels{iPlot}.depths,'LineStyle',...
        allModels{iPlot}.lineStyle,'Color',allModels{iPlot}.color,'DisplayName',...
        allModels{iPlot}.displayName,'LineWidth',line_widths{iPlot});
end
ylim = get(gca,'YLim');
set(gca,'YLim',[ylim(1) ylim(2)]);
set(gca,'YDir','reverse','FontSize',10,'Box','on','Layer','Top');
ax = gca;
xticks(10.^(floor(log10(min(xdata))):2:ceil(log10(max(xdata)))));
xlabel('Resistivity (\Omega-m)');
text(0.90,0.95,char(64+panel_number),'units','normalized','FontSize',14);

end
