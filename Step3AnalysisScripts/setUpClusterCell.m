function KMModels = setUpClusterCell(trueModel,partition,zVals,data,forwardModel)

KMModels = cell(1,partition.numClusters+1);
KMModels{1} = trueModel;
for i = 1:partition.numClusters
    KMModels{i+1} = genModelCalc(10.^partition.centroids(i,:)',zVals,...
        data,rand(1,3),'--',strcat('Centroid #',num2str(i)),forwardModel);
end
end
