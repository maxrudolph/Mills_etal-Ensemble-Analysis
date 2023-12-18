function filenameOut = ensembleAnalysisMaster_noClustering(filename,exact_known)
%{
7/19/21 Ensemble Analysis (no figures)
This represents "step 3" in the process. It takes an ensemble generated in
step 2 and applies a number of analyses, which are then saved in a file for
plotting or further analysis later.
Arguments
filename - a filename for an ensemble file as generated by the inversion script.
exact_known - true if the exact solution is known, otherwise false.
%}

%% Preliminary
addpath(genpath(fileparts(mfilename('fullpath'))))
%access all necessary other folders for scripts
rng(1); %reproducibility
disp('Loading data...')
load(filename,'data','forwardModel','results','pBounds')
slashpos = find(filename == '/',1,'last');
if isempty(slashpos)
    slashpos=1;
end

filenameOut = filename(slashpos+9:end);
disp(['output will be saved to: ' filenameOut]);

%% Model space plots
%The 'model space' plots should show the posterior distribution in
%parameter space (y-axis depth, x-axis rho), as a sort of 3D histogram,
%with any 'appraisals' laid on top.
disp('Model space...')
nzplot = 1000; %number of imaginary (depth)layers to divide appraisals into
nRhoBins = 1001; %number of resistivity bins in model space histogram

%Setup logDepthPlot and logRhoPlot
numSavedRuns = size(results.ensembleRhos,2);
minDistL = -1;%log10(min(data.x));
maxDistL = 5;%log10(max(data.x));
zVals = logspace(minDistL,maxDistL,nzplot)'; %depth values for evaluating
logDepthPlot = log10(repmat(zVals,1,numSavedRuns));
logRhoPlot = zeros(length(zVals),numSavedRuns);
% will contain the (log) depths resistivities
%of every ensemble member, formatted to all be uniform.
residuals = zeros(length(data.x),numSavedRuns);
ewre2n = zeros(1,length(results.ensembleMisfits));
Cdi = inv(data.Cd);
denominator = sqrt(data.y'*Cdi*data.y);
for i = 1:numSavedRuns %for each sln...
    logRhoPlot(:,i) = log10(longForm(zVals,results.ensembleDepths(:,i),...
        results.ensembleRhos(:,i)));
    residuals(:,i) = data.y - forwardModel(results.ensembleDepths(:,i),...
        results.ensembleRhos(:,i),data.lambda);
    ewre2n(i) = sqrt(residuals(:,i)'*Cdi*residuals(:,i))/denominator;
end


%Compute a bivariate histogram of depths/resitvity values which will
%represent the posterior distribution in model space. This will be used not
%only to calculate the maxLikelihood model, but also to plot the
%distribution in model space.
[numElements,binCenters]=hist3([logRhoPlot(:),logDepthPlot(:)],...
    {linspace(log10(pBounds.rhoMin),log10(pBounds.rhoMax),nRhoBins) ...
    log10(zVals)},'CDataMode','auto');
% First linspace is for log(rho), second is for log(depth)

%% Data Space Plots
disp('Data Space...')
nxplot = 500; %Just for plotting. Number of x-values (electrode spacings)
%at which a given solution will be evaluated (put through forward model) at
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

%% Single-model appraisals
disp('Appraisals...')

%Preliminary
meanColor = 'b'; medianColor = 'g'; trueColor = 'r';
bestFitColor = '#df4ec8'; maxLikelihoodColor = '#d1b26f';
msLineStyle = '-'; dsLineStyle = '--';
%These can be changed later

%Model space mean and median
msMeanRhos = 10.^(mean(logRhoPlot,2));%Calculated in log space
mMean = genModelCalc(msMeanRhos,zVals,data,meanColor,msLineStyle,...
    'MS Mean',forwardModel);
msMedianRhos = 10.^(median(logRhoPlot,2));
mMedian = genModelCalc(msMedianRhos,zVals,data,medianColor,msLineStyle,...
    'MS Median',forwardModel);

% Data space median and best fit models + maximumLikelihood
[~,ind2] = sort(ewre2n); %sort all slns by misfit
medianIndex= ind2(floor(length(ind2)/2)); %the one with median misfit
bestIndex = ind2(1);                      %the one with lowest misfit


%Find their corresponding slns in the ensemble
dMedian = genModelInd(medianIndex,zVals,data,medianColor,dsLineStyle,...
    'DS Median',forwardModel,results);
bestFit = genModelInd(bestIndex,zVals,data,bestFitColor,dsLineStyle,...
    'DS Best Fit',forwardModel,results);
% maximumLikelihood = genModelInd()

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

if exact_known
    %Setup true model or exact solution
    [trueDepths,trueRhos] = subStructGen(data.subStructChoice);
    trueDepthsPlot = 10.^logDepthPlot(:,1);
    trueRhoPlot = longForm(trueDepthsPlot,trueDepths,trueRhos);
    trueModel = calculatedModel(trueDepthsPlot,trueRhoPlot,forwardModel(trueDepths,...
        trueRhos,data.lambda),data.y,trueColor,'-','Exact solution');
    trueModel.setWRE2N(data);
else
    trueModel = [];    
end
allModels = {trueModel,mMean,mMedian,maxLikelihood,bestFit,dMedian};

% %% Clustering
% maxNumClusters = 5;
% %You don't have to pick how many clusters you want, the script will find
% %the 'best' number on its own, this is just the maximum you want that
% %number to be. Warning that runtime does not scale linearly with
% %maxNumClusters, so increase at your own risk
% 
% disp('Clustering Euclidean...')
% euclidPartition = clusterMSpace(logRhoPlot,maxNumClusters,'sqeuclidean');
% disp('Clustering Manhattan...')
% manPartition = clusterMSpace(logRhoPlot,maxNumClusters,'cityblock');
% 
% KMModelsEuclid = setUpClusterCell(trueModel,euclidPartition,zVals,...
%     data,forwardModel);
% KMModelsMan = setUpClusterCell(trueModel,manPartition,zVals,...
%     data,forwardModel);
% 
% disp('k-Medoids Euclidean');
% kMedoidsEuclideanPartition = cluster_model_space_kmedoids(logRhoPlot,maxNumClusters,'sqeuclidean');
% disp('k-Medoids Manhattan');
% kMedoidsManhattanPartition = cluster_model_space_kmedoids(logRhoPlot,maxNumClusters,'cityblock');
% 
% disp('Calculating Clustering Models...')
% KMedoidsEuclid    = setUpClusterCell(trueModel,kMedoidsEuclideanPartition,zVals,...
%     data,forwardModel);
% KMedoidsManhattan = setUpClusterCell(trueModel,kMedoidsManhattanPartition,zVals,...
%     data,forwardModel);
% 
% allClusterSets = {KMModelsEuclid, KMModelsMan, KMedoidsEuclid, KMedoidsManhattan};
% allPartitions = {euclidPartition, manPartition, kMedoidsEuclideanPartition, kMedoidsManhattanPartition};
%% Saving
disp('Saving...')

save(['Analysis_' filenameOut],'allModels','binCenters','logRhoPlot',...
    'numElements','nRhoBins','nzplot','xVals','yVals',...
    'zVals','residuals','ewre2n','-v7.3');
%-v7.3 allows for saving of large files
end
