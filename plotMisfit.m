function plotMisfit(inModels,results)
histogram(results.ensembleMisfits,100,'EdgeAlpha',0);
hold on;
yy=get(gca,'YLim');
for iPlot = 1:size(inModels,2)
    plot(inModels{iPlot}.misfit*[1 1],yy,'LineStyle',...
        inModels{iPlot}.lineStyle,'Color',inModels{iPlot}.color,...
        'LineWidth',1.5);
end
set(gca,'FontSize',12);
xlabel('Misfit (\Omega-m)');
<<<<<<< HEAD
=======
text(0.8,0.8,'A','units','normalized','FontSize',14);
>>>>>>> 0ce9c3f5235649a6d740600aa2fbd6288858aff5
end