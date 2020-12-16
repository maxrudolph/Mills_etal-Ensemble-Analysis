function outModel = makeCalculatedModel(inDepths,inRhos,data,...
    forwardModel,colorChoice,title)

    outModel = calculatedModel(inDepths,inRhos,colorChoice,title);
    outModel.setY(forwardModel(outModel.depths,outModel.rhos,...
        data.lambda),data.y);
end
