function ensembleAnalysis3(filename,saveFigures)
rng(1); %reproducibility
disp('Loading data...')
load(filename,'data','forwardModel','results','measure','pBounds')

%% Section 0: Parameter setup
%IF running to generate/save figures, set to true
% saveFigures = true;
nxplot=500; %number of measurement points for evaluating ensemble members
nSavedPlot = 2000; %Number of saved runs to plot
nzplot = 500; %number of imaginary layers to divide models into

meanColor = 'b'; medianColor = 'g';

if saveFigures
    visibility = 'off';
    ensembleName = filename(10:end-9); %captures most relevant info
    slashpos = find(filename == '/',1,'last');
    ensembleName = filename(slashpos+10:end-9);
    folderName = ['figures_' ensembleName];
    mkdir(folderName);
else
    folderName = ' ';
    visibility = 'on';
end

%% Section 1 NF: evaluate ensemble solutions on regularly spaced grid
disp('Evaluating ensemble...');
minDistL = log10(measure.minDist);
maxDistL = log10(measure.maxDist);
numSavedRuns = size(results.ensembleRhos,2);
if nSavedPlot > numSavedRuns
    %If low # of saved runs, plot all, otherwise...
    nSavedPlot = numSavedRuns;
    runPlotIndex = 1:nSavedPlot;
else %... only plot a random subset of them
    runPlotIndex = randperm(numSavedRuns,nSavedPlot);
end
yVals = zeros(nxplot,nSavedPlot); %apparent rhos from measuring at xVals
xVals = logspace(minDistL,maxDistL,nxplot)'; %surface measurement points
zVals = logspace(minDistL,maxDistL,nzplot)'; %depth values for evaluating
lambdaForXVals = makeLambda(xVals); %lambda matrix for forward model
%Ensemble members are saved as media (layers+resistivities), we want to
%show them in data space
for i=1:nSavedPlot
    yVals(:,i) = forwardModel(squeeze(...
        results.ensembleDepths(:,runPlotIndex(i))),...
        results.ensembleRhos(:,runPlotIndex(i)),lambdaForXVals);
end
%yVals is the data space representations of all ensemble members

%% Section 2 NF: Calculating other models
disp('Calculating models...');
logDepthPlot = log10(repmat(zVals,1,numSavedRuns));
logRhoPlot = zeros(size(logDepthPlot));
for i = 1:numSavedRuns %for each run...
    nLayer = nnz(~isnan(results.ensembleDepths(:,i)));
    %...Find the number of layers in that run...
    for j = 1:nLayer %for each layer...
        mask = zVals >= results.ensembleDepths(j,i);
        logRhoPlot(mask,i) = log10(results.ensembleRhos(j,i));
        %Make an appropriate # of values = to that layers resistivity
    end
end
%Generate Model-Space Mean and Median models
inRhos = 10.^(mean(logRhoPlot,2)); %Calculated in log space
[shortDepths,shortRhos] = shortForm(zVals,inRhos);
mMean = calculatedModel(zVals,inRhos,forwardModel(shortDepths,shortRhos,...
    data.lambda),data.y,meanColor,'-','MS Mean');

inRhos = 10.^(median(logRhoPlot,2));
[shortDepths,shortRhos] = shortForm(zVals,inRhos);
mMedian = calculatedModel(zVals,inRhos,forwardModel(shortDepths,shortRhos,...
    data.lambda),data.y,medianColor,'-','MS Median');

% Data space median and best fit models
[~,ind2] = sort(results.ensembleMisfits);
medianIndex= ind2(floor(length(ind2)/2));
bestIndex = ind2(1);
medianRho = results.ensembleRhos(:,medianIndex);
medianDepth = results.ensembleDepths(:,medianIndex);
bestRho = results.ensembleRhos(:,bestIndex);
bestDepth = results.ensembleDepths(:,bestIndex);
medianRhoPlot = zeros(size(zVals,1),1);
bestRhoPlot = zeros(size(zVals,1),1);
for j = 1:size(results.ensembleDepths,1)
    mask1 = zVals >= results.ensembleDepths(j,medianIndex);
    medianRhoPlot(mask1) = results.ensembleRhos(j,medianIndex);
    mask2 = zVals >= results.ensembleDepths(j,bestIndex);
    bestRhoPlot(mask2) = results.ensembleRhos(j,bestIndex);
end
dMedian = calculatedModel(zVals,medianRhoPlot,forwardModel(medianDepth,...
    medianRho,data.lambda),data.y,medianColor,'--','DS Median');
bestFit = calculatedModel(zVals,bestRhoPlot,forwardModel(bestDepth,...
    bestRho,data.lambda),data.y,'#df4ec8','--','DS Best Fit');

%Maximum likelihood model (Model Space)
% compute a bivariate histogram of resitvity values from the
%posterior ensemble
numBins = 2*nxplot; %number of resistivity values for histogram
[numElements,binCenters]=hist3([logRhoPlot(:),logDepthPlot(:)],...
    {linspace(-10,10,numBins) log10(zVals)},'CDataMode','auto');
% First linspace is for log(rho), second is for log(depth)
% at each depth, find the most likely solution (ml_rho)
maxLikelihoodRho = zeros(nzplot,1);
ksRho = linspace(log10(pBounds.rhoMin),log10(pBounds.rhoMax),1e4);
%Bandwidth issues - possible bug in matlab
parfor i=1:nzplot
    % Use ksdensity to approximate the pdf of resistivity at this depth:
    [pdfYVals,pdfXVals] = ksdensity(logRhoPlot(i,:),ksRho,'bandwidth',.05);
    [~,ind1] = max(pdfYVals);
    maxLikelihoodRho(i) = 10.^pdfXVals(ind1);
end
inRhos = maxLikelihoodRho;
[shortDepths,shortRhos] = shortForm(zVals,inRhos);
maxLikelihood = calculatedModel(zVals,inRhos,forwardModel(shortDepths,...
    shortRhos,data.lambda),data.y,'#d1b26f','-','MS Max Likelihood');

%Setup true model/solution
[trueDepths,trueRhos] = modelGen(measure.kMax,measure.modelChoice);
trueDepthsPlot = 10.^logDepthPlot(:,1);
trueRhoPlot = zeros(nzplot,1);
trueNumLayers = nnz(~isnan(trueDepths));
for j = 1:trueNumLayers
    mask = zVals >= trueDepths(j);
    trueRhoPlot(mask) = trueRhos(j);
end
trueModel = calculatedModel(trueDepthsPlot,trueRhoPlot,forwardModel(...
    trueDepths,trueRhos,data.lambda),data.y,'r','-','Exact solution');

%% 3 Figures
disp('Plotting properties...');
smallPlots(results,saveFigures,folderName,visibility);
allModels = {trueModel,mMean,mMedian,maxLikelihood,bestFit,dMedian};
bigPlot(binCenters,numElements,allModels,xVals,yVals,data,results,' ',...
    visibility);
saveFigs(saveFigures,folderName,'4');

%% 4 Clustering stuff
disp('Calculating number of clusters')
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
maxNumClusters = 5;
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
end

function [outDepths,outRhos] = shortForm(inDepths,inRhos)
%removes duplicate values
h = diff(inRhos);
ind = find(h==0);
outRhos = inRhos;
outDepths = inDepths;
outRhos(ind+1) = [];
outDepths(ind+1) = [];
end