function bigPlot(bC,nE,inModels,x,y,data,results,inTitle,saveFigures,...
    folderName,figNumber)
%{
Convenience for making a three-paned plot showing generated models
alongside resistivity distribution (panel A), ensemble members in data
space (panel B), and ensemble members in 
%}
figure('visible',~saveFigures,'units','inches','position',[0 0 6.5 4.5]);
title(inTitle);
subplot(4,2,[2 4 6 8]);
modelSpacePlot(bC,nE,inModels);
%Part 2: Data space plot
subplot(4,2,[3 5 7]);
dataSpacePlot(x,y,inModels,data);

% Part 3 Sub-figure: Histogram of misfit in data space
subplot(4,2,1);
plotMisfit(inModels,results);
saveFigs(saveFigures,folderName,figNumber)

end