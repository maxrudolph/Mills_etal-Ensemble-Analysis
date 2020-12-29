function outModel = makeCalculatedModel(inDepths,inRhos,data,...
    forwardModel,colorChoice,lineStyle,title)

    outModel = calculatedModel(inDepths,inRhos,colorChoice,lineStyle,title);
    outModel.setY(forwardModel(outModel.depths,outModel.rhos,...
        data.lambda),data.y);
end
