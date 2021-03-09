function modelSpacePlot(bC,nE,inModels)
p = pcolor(10.^bC{1},10.^bC{2},nE'); shading flat;
colormap(flipud(bone))
set(gca,'XScale','log','YScale','log','ColorScale','log');
hold on
for i = 1:size(inModels,2)
    plot(inModels{i}.rhos,inModels{i}.depths,'LineStyle',...
        inModels{i}.lineStyle,'Color',inModels{i}.color,'DisplayName',...
        inModels{i}.displayName,'LineWidth',1.75);
end
c=colorbar();
c.Label.String = 'Number of solutions';
set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%legend();
set(gca,'YDir','reverse','FontSize',12,'Box','on','Layer','Top');
k = find(sum(nE')>0);
xlim([min(inModels{1}.rhos)/100,max(inModels{1}.rhos)*100]);
xlabel('Resistivity (\Omega-m)'); ylabel('Depth (m)'); title('Model Space');
lgd = legend('location','northwest');
lgd.FontSize = 8;
text(0.8,0.95,'C','units','normalized','FontSize',14)
end