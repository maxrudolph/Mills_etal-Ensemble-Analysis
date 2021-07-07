function ensembleAnalysisB(logRhoPlot)

disp('Calculating number of clusters')

addpath(genpath(fileparts(mfilename('fullpath'))))
visibility = true; %set false if you don't want figures to appear (as in saving)
rng(1); %reproducibility
disp('Loading data...')
%load(filename,'data','forwardModel','results','pBounds')
numSavedRuns = size(logRhoPlot,2);
maxNumClusters = 5;

%downsample
downsample = false;
if downsample
    downsampleNumber = floor(numSavedRuns/10);
    gmPlotIndex = randperm(numSavedRuns,downsampleNumber);
else
    downsampleNumber = numSavedRuns;
    gmPlotIndex = 1:numSavedRuns;
end

%Use k-means clustering
clust = zeros(downsampleNumber,maxNumClusters);
dsLogRhoPlot = logRhoPlot(:,gmPlotIndex);
stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);
parfor i = 1:maxNumClusters
    fprintf('\nCalculating for %d clusters',i)
    clust(:,i) = kmeans(dsLogRhoPlot',i,'Options',options,'Replicates',5,...
        'MaxIter',1000);
end

%Step 2: Use results from k-means to determine optimal # of clusters
eva = evalclusters(dsLogRhoPlot',clust,'CalinskiHarabasz');
numClusters = eva.OptimalK;

%Step 3: Find Gaussian mixture models
disp('Finding Gaussian models')
GMModel = fitgmdist(dsLogRhoPlot',numClusters,...
    'RegularizationValue',1e-8);

%Step 4: Find kmeans centroids
disp('Finding k-means')
[idxEuclid,CEuclid,sumdEuclid] = kmeans(dsLogRhoPlot',numClusters,...
    'MaxIter',1000,'Replicates',5);
[idxMan,CMan,sumdMan] = kmeans(dsLogRhoPlot',numClusters,'MaxIter',1000,...
    'Distance','cityblock','Replicates',5);

%Step 5: Hierarchical
%disp('Hierarchical clustering')
%T = clusterdata(dsLogRhoPlot','MaxClust',numClusters);

%Step 5: Make models
disp('Calculating models')
GMData = cell(1,numClusters+1);
KMDataEuclid = cell(1,numClusters+1);
KMDataMan = cell(1,numClusters+1);
GMData{1} = trueModel;
KMDataEuclid{1} = trueModel;
KMDataMan{1} = trueModel;
for i = 1:numClusters
    inRhos = 10.^GMModel.mu(i,:)';
    GMData{i+1} = calculatedModel(zVals,inRhos,forwardModel(zVals,inRhos,...
        data.lambda),data.y,rand(1,3),'--',strcat('GM mean ',num2str(i)));
    inRhos = 10.^CEuclid(i,:)';
    KMDataEuclid{i+1} = calculatedModel(zVals,inRhos,forwardModel(zVals,...
        inRhos,data.lambda),data.y,rand(1,3),'--',strcat('Centroid #',...
        num2str(i)));
    inRhos = 10.^CMan(i,:)';
    KMDataMan{i+1} = calculatedModel(zVals,inRhos,forwardModel(zVals,...
        inRhos,data.lambda),data.y,rand(1,3),'--',strcat('Cent ',...
        num2str(i)));
end

%% 5 Plots of GM and k-means
disp('Cluster plotting')

bigPlot(binCenters,numElements,GMData,xVals,yVals,data,results,...
    'Gaussian Mixture Models',visibility)
saveFigs(saveFigures,folderName,'5');
bigPlot(binCenters,numElements,KMDataEuclid,xVals,yVals,data,results,...
    'K-means: Euclidean',visibility);
saveFigs(saveFigures,folderName,'6');
kMeansPlots('K-means: Euclidean',idxEuclid,sumdEuclid,KMDataEuclid,...
    visibility);
saveFigs(saveFigures,folderName,'7');
bigPlot(binCenters,numElements,KMDataMan,xVals,yVals,data,results,...
    'K-means: Manhattan',visibility);
saveFigs(saveFigures,folderName,'8');
kMeansPlots('K-means: Manhattan',idxMan,sumdMan,KMDataMan,visibility);
saveFigs(saveFigures,folderName,'9');

%% 6
if saveFigures
    stringPart1 = 'Ensemble: %s\nFigures recorded on:%s\nTotal runs: %d';%ensembleName,date,numSavedRuns
    stringPart2 = '\n\nMISFITS\nBestFit: %f\nDS Median: %f\nMS Max Likelihood: %f\nMS Mean: %f\nMS Median: %f\nExact solution: %f\n'; %bestFit.misfit,dMedian.misfit,maxLikelihood.misfit,mMean.misfit,mMedian.misfit,trueModel.misfit
    stringPart3 = '\n\nOptimal number of clusters: %d\nMode of ensembleMisfits: %f\n'; %numClusters
    h = histogram(results.ensembleMisfits);
    [~,ind] = max(h.Values);
    readMe = fopen([folderName, '/info.txt'],'w');
    fprintf(readMe,[stringPart1,stringPart2,stringPart3],...
        ensembleName,date,numSavedRuns,bestFit.misfit,...
        dMedian.misfit,maxLikelihood.misfit,mMean.misfit,mMedian.misfit,...
        trueModel.misfit,numClusters,h.BinEdges(ind));
    fclose(readMe);
end
disp('Done');