function modelSpacePanel(binCenters,numElements,allModels,letter)

xdata = 10.^binCenters{1};
mask = xdata >= min(allModels{1}.rhos)/100 & xdata <= max(allModels{1}.rhos)*100;

p = pcolor(xdata(mask),10.^binCenters{2},((1/max(max(numElements)))*(numElements(mask,:)')));
shading flat;

cmap = colormap(flipud(bone));
set(gca,'XScale','log','YScale','log','ColorScale','log');
hold on
for iPlot2 = 1:size(allModels,2)
    plot(allModels{iPlot2}.rhos,allModels{iPlot2}.depths,'LineStyle',...
        allModels{iPlot2}.lineStyle,'Color',allModels{iPlot2}.color,'DisplayName',...
        allModels{iPlot2}.displayName,'LineWidth',1.0);
end
ylim = get(gca,'YLim');
set(gca,'YLim',[ylim(1) ylim(2)]);
set(gca,'YDir','reverse','FontSize',10,'Box','on','Layer','Top');
ax = gca;
xticks(10.^(floor(log10(min(xdata))):2:ceil(log10(max(xdata)))));
xlabel('Resistivity (\Omega-m)');
text(0.05,0.95,letter,'units','normalized','FontSize',14);

end
