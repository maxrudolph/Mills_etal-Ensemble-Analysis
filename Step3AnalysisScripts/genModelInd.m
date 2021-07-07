function outModel = genModelInd(index,zVals,data,color,lineStyle,title,...
    forwardModel,results)

rhos = results.ensembleRhos(:,index);
depths = results.ensembleDepths(:,index);
rhoPlot = 10.^longForm(zVals,depths,rhos);
outModel = calculatedModel(zVals,rhoPlot,...
    forwardModel(depths,rhos,data.lambda),data.y,color,lineStyle,title);
end