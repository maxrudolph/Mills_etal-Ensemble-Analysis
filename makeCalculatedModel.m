function outModel = makeCalculatedModel(inDepths,inRhos,data,...
    forwardModel,colorChoice)

    outModel = calculatedModel(inDepths,inRhos,colorChoice);
    outModel.setY(forwardModel(outModel.depths,outModel.rhos,...
        data.lambda),data.y);
end
