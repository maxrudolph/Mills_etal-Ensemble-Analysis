function h = modelSpacePanel(binCenters,numElements,allModels,panel_number,line_widths,piecewiseLinear,pBounds)

xdata = 10.^binCenters{1};
if ~isempty(allModels) && isfield(allModels{1},'rhos')
    mask = xdata >= min(allModels{1}.rhos)/100 & xdata <= max(allModels{1}.rhos)*100;
else
    % mask = xdata >= 5e-1 & xdata <= 5e4;
    mask = true(size(xdata));
end
% Pseudocolor plot of model-space probability density
p = pcolor(xdata(mask),10.^binCenters{2},((1/max(max(numElements)))*(numElements(mask,:)')));
shading flat;
cmap = colormap(flipud(bone));
set(gca,'XScale','log','YScale','log','ColorScale','log');
hold on

% add curves for model space solutions
for iPlot = 1:length(allModels) % iterate backwards to plot exact solution last.
    if (isfield(allModels{iPlot},'rhos')) ...
            || (isobject(allModels{iPlot}) && isprop(allModels{iPlot},'rhos'))
        % interpolate solution if this is a piecewiseLinear model
        if piecewiseLinear && iPlot>1 % don't interpolate the exact solution.
            [depth_plot,rho_plot] = piecewiseLinearSolution(allModels{iPlot}.depths,allModels{iPlot}.rhos,pBounds);
        else
            depth_plot = allModels{iPlot}.depths;
            rho_plot = allModels{iPlot}.rhos;
        end
        plot(rho_plot,depth_plot,'LineStyle',...
            allModels{iPlot}.lineStyle,'Color',allModels{iPlot}.color,'DisplayName',...
            allModels{iPlot}.displayName,'LineWidth',line_widths{iPlot});
    end
end
ylim = get(gca,'YLim');
set(gca,'YLim',[ylim(1) ylim(2)]);
set(gca,'YDir','reverse','FontSize',10,'Box','on','Layer','Top');
ax = gca;
xticks(10.^(floor(log10(min(xdata))):2:ceil(log10(max(xdata)))));
xlabel('Resistivity (\Omega-m)');
text(0.90,0.95,char(64+panel_number),'units','normalized','FontSize',14);

end
