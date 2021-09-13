function KMModels = setUpClusterCell(trueModel,partition,zVals,data,forwardModel)


C = [0 0 0;
    colororder()
    ];
KMModels = cell(1,partition.numClusters+1);
KMModels{1} = trueModel;
for i = 1:partition.numClusters
    KMModels{i+1} = genModelCalc(10.^partition.centroids(i,:)',zVals,...
        data,C(i,:),'--',strcat('Centroid #',num2str(i)),forwardModel);
end
end
