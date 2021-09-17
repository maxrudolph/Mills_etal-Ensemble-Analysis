function outputPartition = cluster_model_space_kmedoids(lRPlot,maxNum,distMetric)
% Compute cluster analysis in model space.
% lRPlot is the logarithm of resistivity and has dimensions:
%      (number of depths) x (number of solutions).
% maxNum is the maximum number of clusters to be computed
% distMetric is the distance method to use.
idxs = zeros(size(lRPlot,2),maxNum);
cents = zeros(sum(1:maxNum),size(lRPlot,1));
sumds = zeros(sum(1:maxNum),1);

stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream,'MaxIter',1000,'Display','iter');
j = 1;
for i = 1:maxNum
    fprintf('Calculating for %d clusters\n',i)
    index = j:(j+i-1);
    [idxs(:,i),cents(index,:),sumds(index)] = kmedoids(lRPlot',i,...
        'Options',options,'Replicates',10,'Distance',distMetric);
    j = j+i;
end

eva = evalclusters(lRPlot',idxs,'CalinskiHarabasz');
numClusters = eva.OptimalK;

index = 1+sum(1:numClusters-1);
outputPartition.numClusters = numClusters;
outputPartition.indices = idxs(:,numClusters);
outputPartition.centroids = cents(index:index+numClusters-1,:);
outputPartition.sumd = sumds(index:index+numClusters-1);
end
