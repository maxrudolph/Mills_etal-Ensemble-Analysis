function ensembleAnalysisA(filename)
%{
6 /23/2021 Ensemble Analysis, the last step
This script starts the process by setting up parameters, evaluating slns,
and calculating single-model appraisals.
Inputs: filename linking to a file produced by the 'inversion' script. 



%}
rng(1); %reproducibility
disp('Loading data...')
load(filename,'data','forwardModel','results','pBounds')

%% Section 0: Parameter setup
%IF running to generate/save figures, set to true
% saveFigures = true;
nxplot=500; %number of measurement points for evaluating ensemble members
nSavedPlot = 2000; %Number of saved runs to plot (not evaluate, just plot)
nzplot = 500; %number of imaginary layers to divide appraisals into
%ie, resolution

meanColor = 'b'; medianColor = 'g';



%% Section 1 NF: evaluate ensemble solutions on regularly spaced grid
disp('Evaluating ensemble...');
minDistL = log10(min(data.x));
maxDistL = log10(max(data.x));
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
lambdaForXVals = makeLambda(xVals,size(data.lambda,1)); 
%lambda matrix for forward model, second parameter is filter size
%Ensemble members are saved as subsurface structures (layers+resistivities), 
%we want to show them in data space
for i=1:nSavedPlot
    yVals(:,i) = forwardModel(squeeze(...
        results.ensembleDepths(:,runPlotIndex(i))),...
        results.ensembleRhos(:,runPlotIndex(i)),lambdaForXVals);
end
%yVals is the data space representations of selected ensemble members

%% Section 2: Format ensemble solutions for evaluation
disp('Calculating models...');
logDepthPlot = log10(repmat(zVals,1,numSavedRuns));
logRhoPlot = zeros(size(logDepthPlot));
%logDepthPlot and logRhoPlot will contain the depths and resistivities of
%every ensemble member, formatted to all be uniform.
for i = 1:numSavedRuns %for each run...
    nLayer = nnz(~isnan(results.ensembleDepths(:,i)));
    %...Find the number of layers in that run...
    for j = 1:nLayer %for each layer...
        mask = zVals >= results.ensembleDepths(j,i);
        logRhoPlot(mask,i) = log10(results.ensembleRhos(j,i));
        %Make an appropriate # of values = to that layers resistivity
    end
end

%% Section 3: Calculate models
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

function [outDepths,outRhos] = shortForm(inDepths,inRhos)
%removes duplicate values
h = diff(inRhos);
ind = find(h==0);
outRhos = inRhos;
outDepths = inDepths;
outRhos(ind+1) = [];
outDepths(ind+1) = [];
end