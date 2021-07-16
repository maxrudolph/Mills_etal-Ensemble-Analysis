function ensembleAnalysisMaster(filename, saveFigures)

%% Preliminary
addpath(genpath(fileparts(mfilename('fullpath'))))
%access all necessary other folders for scripts
rng(1); %reproducibility
disp('Loading data...')
load(filename,'data','forwardModel','results','pBounds')

if saveFigures
    slashpos = find(filename == '/',1,'last');
    ensembleName = filename(slashpos+10:end-9);
    folderName = ['figures_' ensembleName];
    mkdir(folderName);
else
    folderName = ' ';
end

%% Figures 1 through 3
disp('Producing first three figures...')
firstThreeFigures(results,saveFigures,folderName)

%% Model space plots (figures 4,5,7)
disp('Model space...')
nzplot = 1000; %number of imaginary (depth)layers to divide appraisals into
nRhoBins = 1000; %number of resistivity bins in model space histogram

%Setup logDepthPlot and logRhoPlot
numSavedRuns = size(results.ensembleRhos,2);
minDistL = log10(min(data.x));
maxDistL = log10(max(data.x));
zVals = logspace(minDistL,maxDistL,nzplot)'; %depth values for evaluating
logDepthPlot = log10(repmat(zVals,1,numSavedRuns));
logRhoPlot = zeros(length(zVals),numSavedRuns);
% will contain the (log) depths resistivities
%of every ensemble member, formatted to all be uniform.
for i = 1:numSavedRuns %for each sln...
    logRhoPlot(:,i) = log10(longForm(zVals,results.ensembleDepths(:,i),...
        results.ensembleRhos(:,i)));
end

%Compute a bivariate histogram of depths/resitvity values which will
%represent the posterior distribution in model space. This will be used not
%only to calculate the maxLikelihood model, but also to plot the
%distribution in model space.
[numElements,binCenters]=hist3([logRhoPlot(:),logDepthPlot(:)],...
    {linspace(log10(pBounds.rhoMin),log10(pBounds.rhoMax),nRhoBins) ...
    log10(zVals)},'CDataMode','auto');
% First linspace is for log(rho), second is for log(depth)

%% Data Space Plots for figs 4,5,7
disp('Data Space...')
nxplot = 500;
filterSize = size(data.lambda,1);
nSavedPlot = 2000; %number of runs to plot
numSavedRuns = size(results.ensembleRhos,2);

if nSavedPlot >= numSavedRuns
    %If low # of saved runs, plot all, otherwise...
    nSavedPlot = numSavedRuns;
    runPlotIndex = 1:nSavedPlot;
else %... only plot a random subset of them
    runPlotIndex = randperm(numSavedRuns,nSavedPlot);
end
yVals = zeros(nxplot,nSavedPlot); %apparent rhos from measuring at xVals
xVals = logspace(minDistL,maxDistL,nxplot)'; %surface measurement points
lambdaForXVals = makeLambda(xVals,filterSize);
%lambda matrix for forward model
%Ensemble members are saved as subsurface structures (layers+resistivities),
%we want to show them in data space
for i=1:nSavedPlot
    yVals(:,i) = forwardModel(squeeze(...
        results.ensembleDepths(:,runPlotIndex(i))),...
        results.ensembleRhos(:,runPlotIndex(i)),lambdaForXVals);
end
%yVals is the data space representations of selected ensemble members
%remember, only the ones getting plotted, not all ensemble members

%% Single-model appraisals for figure 4
disp('Appraisals...')

%Preliminary
meanColor = 'b'; medianColor = 'g'; trueColor = 'r';
bestFitColor = '#df4ec8'; maxLikelihoodColor = '#d1b26f';
msLineStyle = '-'; dsLineStyle = '--';

%Model space mean and median
msMeanRhos = 10.^(mean(logRhoPlot,2));%Calculated in log space
mMean = genModelCalc(msMeanRhos,zVals,data,meanColor,msLineStyle,...
    'MS Mean',forwardModel);
msMedianRhos = 10.^(median(logRhoPlot,2));
mMedian = genModelCalc(msMedianRhos,zVals,data,medianColor,msLineStyle,...
    'MS Median',forwardModel);

% Data space median and best fit models
[~,ind2] = sort(results.ensembleMisfits); %sort all slns by misfit
medianIndex= ind2(floor(length(ind2)/2)); %the one with median misfit
bestIndex = ind2(1);                      %the one with lowest misfit
%Find their corresponding slns in the ensemble
dMedian = genModelInd(medianIndex,zVals,data,medianColor,dsLineStyle,...
    'DS Median',forwardModel,results);
bestFit = genModelInd(bestIndex,zVals,data,bestFitColor,dsLineStyle,...
    'DS Best Fit',forwardModel,results);

%Maximum likelihood model (Model Space)
% At each depth, find the most likely solution
maxLikelihoodRho = zeros(nzplot,1);
ksRho = linspace(log10(pBounds.rhoMin),log10(pBounds.rhoMax),1e4);
%This loop will take up ~85-90% of runtime if not using parallel computing
parfor i=1:nzplot
    % Use ksdensity to approximate the pdf of resistivity at this depth:
    [pdfYVals,pdfXVals] = ksdensity(logRhoPlot(i,:),ksRho,'bandwidth',.05);
    %Bandwidth issues - possible bug in matlab?
    [~,ind1] = max(pdfYVals); %Find index of highest probability
    maxLikelihoodRho(i) = 10.^pdfXVals(ind1); %Find associated resistivity
end
maxLikelihood = genModelCalc(maxLikelihoodRho,zVals,data,...
    maxLikelihoodColor,msLineStyle,'MS Max Likelihood',forwardModel);

%Setup true model/solution
[trueDepths,trueRhos] = subStructGen(data.subStructChoice);
trueDepthsPlot = 10.^logDepthPlot(:,1);
trueRhoPlot = longForm(trueDepthsPlot,trueDepths,trueRhos);
trueModel = genModelCalc(trueRhoPlot,trueDepthsPlot,data,trueColor,'-',...
    'Exact solution',forwardModel);

allModels = ...
    {trueModel,mMean,mMedian,maxLikelihood,bestFit,dMedian};

%% Figure 4
disp('Figure 4...')
bigPlot(binCenters,numElements,allModels,xVals,yVals,data,results,' ',...
    saveFigures,folderName,'4');


%% Clustering for Figures 5 - 6

maxNumClusters = 8;

disp('Clustering Euclidean...')
euclidPartition = clusterMSpace(logRhoPlot,maxNumClusters,'sqeuclidean');
disp('Clustering Manhattan...')
manPartition = clusterMSpace(logRhoPlot,maxNumClusters,'cityblock');

disp('Calculating Clustering Models...')
KMModelsEuclid = cell(1,euclidPartition.numClusters+1);
KMModelsMan = cell(1,manPartition.numClusters+1);
KMModelsEuclid{1} = trueModel;
KMModelsMan{1} = trueModel;
for i = 1:euclidPartition.numClusters
    KMModelsEuclid{i+1} = genModelCalc(...
        10.^euclidPartition.centroids(i,:)',zVals,data,rand(1,3),'--',...
        strcat('Centroid #',num2str(i)),forwardModel);
end
for i = 1:manPartition.numClusters
    KMModelsMan{i+1} = genModelCalc(...
        10.^manPartition.centroids(i,:)',zVals,data,rand(1,3),'--',...
        strcat('Centroid #',num2str(i)),forwardModel);
end


%% 5 Plots of GM and k-means
disp('Cluster plotting...')

bigPlot(binCenters,numElements,KMModelsEuclid,xVals,yVals,data,results,...
    'K-means: Euclidean',saveFigures,folderName,'5');
bigPlot(binCenters,numElements,KMModelsMan,xVals,yVals,data,results,...
    'K-means: Manhattan',saveFigures,folderName,'6');

disp('Saving...')
filenameOut = ['Analysis_',filename(9:end)];

save(filenameOut,'allModels','binCenters','euclidPartition',...
    'KMModelsEuclid','KMModelsMan','logRhoPlot','manPartition',...
    'nRhoBins','nzplot','zVals','-v7.3'); %-v7.3 allows for saving of large files
end
