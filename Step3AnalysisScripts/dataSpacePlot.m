function dataSpacePlot(x,y,inModels,data)
hold on;
ensembleColor = [200 200 200]/255;

%Ensemble solution
for i=1:size(y,2)-1
    plot(x,y(:,i),'Color',ensembleColor);
end
i = size(y,2);
h = plot(x,y(:,i),'Color',ensembleColor,'DisplayName','Ensemble Members');

set(gcf,'Color','w');
for i = 1:size(inModels,2)
    plot(data.x,inModels{i}.y,'Color',inModels{i}.color,'LineStyle',...
        inModels{i}.lineStyle,'LineWidth',1.0,...
        'DisplayName',inModels{i}.displayName);
end
%h1 = plot(x,mean(y,2),'b--','LineWidth',1,'DisplayName','DS Mean');
h2 = plot(data.x,data.y,'.','Color',inModels{1}.color,'MarkerSize',10.0,...
    'DisplayName','Data + noise');
lgd = legend([h,h2],'Location','best');
lgd.FontSize = 7;
set(gca,'FontSize',10,'Color','w','XScale','log','YScale','log','Box','on');
xlabel('Array Spacing (m)'); ylabel('Apparent Resistivity (\Omega-m)')
%text(0.9,0.95,'B','units','normalized','FontSize',14)
%title('Data Space')
end