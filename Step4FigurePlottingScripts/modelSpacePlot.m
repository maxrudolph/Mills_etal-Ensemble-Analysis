function modelSpacePlot(bC,nE,inModels)
xdata = 10.^bC{1};
mask = xdata >= min(inModels{1}.rhos)/100 & xdata <= max(inModels{1}.rhos)*100;
% xlim([min(inModels{1}.rhos)/100,max(inModels{1}.rhos)*100]);
xticks = 10.^(floor(log10(min(xdata))):2:ceil(log10(max(xdata))));
p = pcolor(xdata(mask),10.^bC{2},nE(mask,:)'); shading flat;

colormap(flipud(bone))
set(gca,'XScale','log','YScale','log','ColorScale','log');
hold on
for i = 1:size(inModels,2)
    plot(inModels{i}.rhos,inModels{i}.depths,'LineStyle',...
        inModels{i}.lineStyle,'Color',inModels{i}.color,'DisplayName',...
        inModels{i}.displayName,'LineWidth',1.0);
end
ylim = get(gca,'YLim');
set(gca,'YLim',[ylim(1) ylim(2)]);
c=colorbar();
c.Label.String = 'Number of solutions';
set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend();
set(gca,'YDir','reverse','FontSize',10,'Box','on','Layer','Top');
k = find(sum(nE')>0);

ax = gca;
ax.XTickMode = 'manual';
ax.XTick = xticks;
xlabel('Resistivity (\Omega-m)'); ylabel('Depth (m)');
title('Solution Space');
lgd = legend('location','southeast');
lgd.FontSize = 7;
%text(0.9,0.95,'C','units','normalized','FontSize',14)
end