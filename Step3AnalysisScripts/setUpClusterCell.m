function KMModels = setUpClusterCell(trueModel,partition,zVals,data,forwardModel)


C = [0 0 0;
    colororder()
    ];

if ~isempty(trueModel)
    KMModels = cell(1,partition.numClusters+1);
    KMModels{1} = trueModel;
    offt = 1;
else
    KMModels = cell(1,partition.numClusters);
    offt = 0;
end

for i = 1:partition.numClusters
    KMModels{i+offt} = genModelCalc(10.^partition.centroids(i,:)',zVals,...
        data,C(i,:),'--',strcat('Centroid #',num2str(i)),forwardModel);
end
end
