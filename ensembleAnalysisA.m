function ensembleAnalysisA(filename,saveFigures)
%{
7/6/2021 Ensemble Analysis, the last step
This script starts the process by setting up parameters, evaluating slns,
and generating single-model appraisals: Model-space mean and median,
data-space median, best fit, and maximum likelihood (plus 'exact sln')
Inputs: filename linking to a file produced by the 'inversion' script.
This will contain structures:
    data: a structure with fields:
        x: a vector of x-values, or electrode spacings
        lambda: a 'lambda matrix' generated from x, see 'makeLambda' script
        fx: a vector of y-values (apparent resistivity measurements) at the
            corresponding x-values from above, with NO noise added
        y: vector of y-values WITH noise
        Cd: covariance matrix
        subStructChoice: string indicating which subsurface structure was
            used to generate the data, see subStructGen script
        noiseCoef: how noisy the data was
    results: a structure with fields
        ensembleDepths, ensembleRhos, ensembleVars, ensembleMisfits,
            ensembleNumLayers: the specifics of every "sln" (subsurface
            structure) in the ensemble. Depths of each layer interface (m),
            resistivity of each layer (ohm-meters), variance, misfit with
            data, and number of layers, respectively
        allChoices: what random choice was made at each iteration of MCMC
        allLikelihoods, allMisfits, allProbAccepts, allVars: the
            likelihood, misfit, acceptance probability and variance of ALL
            solutions, not just those that were saved
        maxLayers: The maximum # of layers allowed at each iter by burn-in
    pBounds: structure containing parameter bounds for everything, but
        really we only need min and max resistivity bounds
    forwardModel: function handle of what the forward model used in the
        inversion was
%}
addpath(genpath(fileparts(mfilename('fullpath'))))
%access all necessary other folders for scripts
rng(1); %reproducibility
disp('Loading data...')
load(filename,'data','forwardModel','results','pBounds')

if saveFigures
    visibility = 'off'; %Figures won't appear
    ensembleName = filename(10:end-9); %captures most relevant info
    slashpos = find(filename == '/',1,'last');
    ensembleName = filename(slashpos+10:end-9);
    folderName = ['figures_' ensembleName];
    mkdir(folderName);
else
    folderName = ' ';
    visibility = 'on';
end

%% Section 0: Parameter setup

nxplot=500; %number of measurement points for evaluating ensemble members
%nxplot is only used for the ensemble members that get plotted
nSavedPlot = 2000; %Number of saved runs to plot (not evaluate, just plot)
nzplot = 1000; %number of imaginary (depth)layers to divide appraisals into
nRhoBins = 1000; %number of resistivity bins in model space histogram

%set colors and line styles for plotting
meanColor = 'b'; medianColor = 'g'; trueColor = 'r';
bestFitColor = '#df4ec8'; maxLikelihoodColor = '#d1b26f';
msLineStyle = '-'; dsLineStyle = '--';

filterSize = size(data.lambda,1);


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
%yVals and xVals are only used for plotting, zVals is used later
zVals = logspace(minDistL,maxDistL,nzplot)'; %depth values for evaluating
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

%% Section 2: Format ensemble solutions for evaluation
disp('Formatting...');
logDepthPlot = log10(repmat(zVals,1,numSavedRuns));
logRhoPlot = generateLogRhoPlot(results.ensembleDepths,...
    results.ensembleRhos,zVals); 
%resistivities at each logDepthPlot depth for each ensemble member

%% Section 3: Calculate models
%Generate Model-Space Mean and Median models
disp('Calculating models...');

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
% compute a bivariate histogram of depths/resitvity values which will 
%represent the posterior distribution in model space. This will be used not
%only to calculate the maxLikelihood model, but also to plot the
%distribution in model space.
[numElements,binCenters]=hist3([logRhoPlot(:),logDepthPlot(:)],...
    {linspace(log10(pBounds.rhoMin),log10(pBounds.rhoMax),nRhoBins) ...
    log10(zVals)},'CDataMode','auto');
% First linspace is for log(rho), second is for log(depth)
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
%{
analyzedEnsemble.binCenters = binCenters;
analyzedEnsemble.numElements = numElements;
analyzedEnsemble.allModels = ...
    {trueModel,mMean,mMedian,maxLikelihood,bestFit,dMedian};
analyzedEnsemble.xVals = xVals;
analyzedEnsemble.yVals = yVals;
%}
    allModels = ...
    {trueModel,mMean,mMedian,maxLikelihood,bestFit,dMedian};
smallPlots(results,visibility,saveFigures,folderName);
bigPlot(binCenters,numElements,allModels,xVals,yVals,data,results,' ',...
    visibility);
saveFigs(saveFigures,folderName,'4');





end