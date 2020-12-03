%function ensembleAnalysis(filename)
disp('Loading data...')
load(filename,'data','forwardModel','results','measure','pBounds')

%% Section 1 NF: evaluate ensemble solutions on regularly spaced grid
disp('Evaluating ensemble...');
minDistL = log10(measure.minDist);
maxDistL = log10(measure.maxDist);
nxplot=200; %number of measurement points
nSavedPlot = 2000; %Number of saved runs to plot
if nSavedPlot > size(results.ensembleRhos,2)
    %If low # of saved runs, plot all, otherwise...
    nSavedPlot = size(results.ensembleRhos,2);
    runPlotIndex = 1:nSavedPlot;
else %... only plot a random subset of them
    runPlotIndex = randperm(size(results.ensembleRhos,2),nSavedPlot);
end
yVals = zeros(nxplot,nSavedPlot);
xVals = logspace(minDistL,maxDistL,nxplot)';
lambdaForXVals = makeLambda(xVals);

for i=1:nSavedPlot
    yVals(:,i) = forwardModel(squeeze(results.ensembleDepths(:,runPlotIndex(i))),...
        results.ensembleRhos(:,runPlotIndex(i)),lambdaForXVals);
end

%% Section 2 NF: Calculating other models
disp('Calculating models...');
%Calculated mean and calculated median
%Both calculated in log space, from the entire ensemble
nzplot = 500;
zVals = 10.^linspace(log10(pBounds.depthMin),log10(pBounds.depthMax),nzplot)';
numSavedRuns = size(results.ensembleRhos,2);
depthPlot = repmat(zVals,1,numSavedRuns);
logDepthPlot = log10(depthPlot);
rhoPlot = zeros(size(depthPlot));
for i = 1:numSavedRuns %for each run...
    nLayer = nnz(~isnan(results.ensembleDepths(:,i)));
    %...Find the number of layers in that run...
    for j = 1:nLayer
        mask = zVals >= results.ensembleDepths(j,i);
        rhoPlot(mask,i) = results.ensembleRhos(j,i);
    end
end
logRhoPlot = log10(rhoPlot);
calcMean = makeCalculatedModel(zVals,10.^(mean(logRhoPlot,2)),data,...
    forwardModel,'g');
calcMedian = makeCalculatedModel(zVals,10.^(median(logRhoPlot,2)),data,...
    forwardModel,'y');

% Ensemble median and best fit models
[~,ind2] = sort(results.ensembleMisfits);
medianIndex= ind2(floor(length(ind2)/2));
bestIndex = ind2(1);
ensembleMedian = makeCalculatedModel(...
    results.ensembleDepths(:,medianIndex),...
    results.ensembleRhos(:,medianIndex),data,forwardModel,'k');
bestFit = makeCalculatedModel(results.ensembleDepths(:,bestIndex),...
    results.ensembleRhos(:,bestIndex),data,forwardModel,'c');

%Maximum likelihood model
% compute a bivariate histogram of resitvity values from the posterior ensemble
numBins = 2*nxplot;
[numElements,binCenters]=hist3([logRhoPlot(:),logDepthPlot(:)],...
    {linspace(-10,10,numBins) log10(zVals)},'CDataMode','auto');
% First linspace is for log(rho), second is for log(depth)
% at each depth, find the most likely solution (ml_rho)
maxLikelihoodRho = zeros(nzplot,1);
ksRho = linspace(log10(pBounds.rhoMin),log10(pBounds.rhoMax),1e4);
%Bandwidth issues - possible bug in matlab
logRhoPlot = logRhoPlot';
parfor i=1:nzplot
    i
    % Use ksdensity to approximate the pdf of resistivity at this depth:
    [pdfYVals,pdfXVals] = ksdensity(logRhoPlot(:,i),ksRho,'bandwidth',.01);
    [~,ind1] = max(pdfYVals);
    maxLikelihoodRho(i) = 10.^pdfXVals(ind1);
end
logRhoPlot = logRhoPlot';
maxLikelihood = makeCalculatedModel(zVals,maxLikelihoodRho,data,...
    forwardModel,'m');

%Setup true model 
[trueDepths,trueRhos] = modelGen(measure.kMax,measure.modelChoice);
trueModel = makeCalculatedModel(trueDepths,trueRhos,data,forwardModel,...
    'r');

%% 3 Figures
disp('Plotting properties...');
smallPlots(results);
%% Section 4 The big figure:
disp('More plots...');

% Part 1: Data space plot
figure;
subplot(4,2,[3 5 7]);
hold on;
ensembleColor = [200 200 200]/255;

%Ensemble solution
for i=1:nSavedPlot
    hEnsemble=plot(xVals,yVals(:,i),'Color',ensembleColor);
end
set(gca,'Box','on');
set(gcf,'Color','w');
allModels = {trueModel,calcMean,calcMedian,bestFit,maxLikelihood,...
    ensembleMedian};
plots = cell(1,6);
for iPlot = 1:6
    plots{iPlot} = loglog(data.x,allModels{iPlot}.y,'-','Color',...
        allModels{iPlot}.color);
end

hdata = loglog(data.x,data.y,'r.','MarkerSize',10.0);
hEnsembleMean = loglog(xVals,mean(yVals,2),'k*');

legend([plots{1},hdata,hEnsemble,hEnsembleMean,plots{6},plots{2},...
    plots{3},plots{4},plots{5}]...
    ,{'Exact','Data+noise','Ensemble','Ensemble mean','Ensemble median',...
    'Calculated mean','Calculated median','Best fit model',...
    'Max Likelihood Model'},'Position',[0.2 0.2 0.1 0.2]);

set(gca,'FontSize',12);
set(gca,'Color','w');
set(gca,'XScale','log','YScale','log');
xlabel('Array Spacing (m)');
ylabel('Apparent Resistivity (\Omega-m)')

% Part 2 Sub-figure: Model space plot
subplot(4,2,[2 4 6 8]);
pcolor(10.^binCenters{1},10.^binCenters{2},numElements'); shading flat;
set(gca,'XScale','log','YScale','log');

hold on
trueLogRhoPlot = zeros(nzplot,1);
trueLogDepthsPlot = logDepthPlot(:,1);

trueNumLayers = nnz(~isnan(trueDepths));
for j=1:trueNumLayers
    mask = log10(zVals) >= log10(trueDepths(j));
    trueLogRhoPlot(mask) = log10(trueRhos(j));
end
plot(10.^trueLogRhoPlot,10.^trueLogDepthsPlot,trueModel.color);
plot(rhoPlot(:,medianIndex),10.^logDepthPlot(:,medianIndex),'-',...
    'Color',ensembleMedian.color);
plot(rhoPlot(:,bestIndex),10.^logDepthPlot(:,bestIndex),'-',...
    'Color',bestFit.color);
plot(maxLikelihood.rhos,maxLikelihood.depths,'-',...
    'Color',maxLikelihood.color);
plot(calcMean.rhos,calcMean.depths,'-','Color',calcMean.color);
colorbar();
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'Box','on');
xlabel('Resistivity (\Omega-m)');
ylabel('Depth (m)');


% Part 3 Sub-figure: Histogram of misfit in data space
subplot(4,2,1);
histogram(results.ensembleMisfits,100);
hold on;
yy=get(gca,'YLim');
plot(calcMedian.misfit*[1 1],yy,'-','Color',calcMedian.color);
plot(bestFit.misfit*[1 1],yy,'-','Color',bestFit.color);
plot(calcMean.misfit*[1 1],yy,'-','Color',calcMean.color);
plot(maxLikelihood.misfit*[1 1],yy,'-','Color',maxLikelihood.color);
plot(ensembleMedian.misfit*[1 1],yy,'*','Color',ensembleMedian.color);
set(gca,'FontSize',12);
xlabel('Misfit (m)');
f=gcf;

%% 5 Optional: Save ensemble for The Sequencer
% nsequence = 100000;
% slashpos = find(filename=='/',1,'last');
% txtfile = [filename(slashpos+1:end-3) 'csv'];
% % each row should correspond to a different resistivity profile
% csvwrite(txtfile,logRhoPlot(:,end-nsequence:end)'); % note that sequencer
% struggles with very large datasets. I choose just the last nsequence
% samples in the ensemble.

%% 6 Optional: save the figure
% export_fig([results_file(1:end-4) '_ensemble.eps'])
%exportgraphics(gcf,[results_file(1:end-4) '_ensemble.eps'],'ContentType','vector');

%end


%% Figure of model-space vs data-space misfit
%{
allModelSpaceMisfits = zeros(numSavedRuns,1);
for i = 1:numSavedRuns
    allModelSpaceMisfits(i) = norm(logRhoPlot(:,i)-...
        log10(interpolatedTrueRhos'))/norm(log10(interpolatedTrueRhos'));
end
allDataSpaceMisfits = results.ensembleMisfits/norm(data.y);
figure,plot(allModelSpaceMisfits,allDataSpaceMisfits,'.')
xlabel('Model Space Misfit')
ylabel('Data Space Misfit')
figure
hold on
[~,ind] = sort(allModelSpaceMisfits);
for i = 1:10
    plot(logRhoPlot(:,ind(end-i+1)))
    plot(logRhoPlot(:,ind(i)),'--')
end

%% Next figure

disp('kmeans');
numClusters = 2;
[idx,C,sumd] = kmeans(logRhoPlot',numClusters);
figure;
pcolor(10.^c{1},10.^c{2},N'); shading flat;
set(gca,'XScale','log','YScale','log');
hold on
for i = 1:size(C,1)
    plot(10.^C(i,:),10.^logDepthPlot(:,1),'--','LineWidth',2);
end
plot(10.^trueLogRhoPlot,10.^trueLogDepthsPlot,trueColor);
colorbar();
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'Box','on');
xlabel('Resistivity (\Omega-m)');
ylabel('Depth (m)');
legend()

figure();
subplot(3,2,1);
hist(idx);
%set(gca,'YScale','log');
title('Number of ensemble slns');
xlabel('Cluster #');
ylabel('Frequency');

numRunsEachCluster = zeros(1,numClusters);
averageDist = zeros(1,numClusters);
centroidMisfits = zeros(1,numClusters);
centroidYs = zeros(length(data.y),numClusters);
for i = 1:numClusters
    numRunsEachCluster(i) = nnz(idx==i);
    averageDist(i) = sumd(i)/numRunsEachCluster(i);
    centroidYs(:,i) = forwardModel(xVals,10.^C(i,:),data.lambda);
    % evaluate the forward model for the maximum likelihood.;
    centroidMisfits(i) = norm(centroidYs(:,i)-data.y);
end

subplot(3,2,2);
bar(averageDist);
title('Distance from Centroid ');
ylabel('Log-distance');
xlabel('Cluster');

subplot(3,2,3);
bar(centroidMisfits)
title('Misfit')
ylabel('Data-space misfit (m)')
xlabel('Cluster');

subplot(3,2,4);
histogram(results.ensembleMisfits,100);
hold on;
yy=get(gca,'YLim');
for i = 1:numClusters
    plot(centroidMisfits(i)*[1 1],yy,'LineWidth',1);
end
set(gca,'FontSize',12);
xlabel('Log misfit');
f=gcf;
plot(medianModelMisfit*[1 1],yy,'Color',medianColor,'LineWidth',1);
plot(bestFitModelMisfit*[1 1],yy,'Color',bestFitColor,'LineWidth',1);
plot(meanModelMisfit*[1 1],yy,'Color',meanColor,'LineWidth',1);
plot(maxLikelihoodMisfit*[1 1],yy,'Color',maxLikelihoodColor,'LineWidth',1);
plot(ensembleMedianModelMisfit*[1 1],yy,ensembleMedianColor);


subplot(3,2,5);
hold on;

%Ensemble solution
for i=1:nSavedPlot
    hEnsemble=plot(xVals,yVals(:,i),'Color',ensembleColor);
end
set(gca,'Box','on');
set(gcf,'Color','w');
%True model
hdata = loglog(data.x,data.y,'r.','MarkerSize',10.0);
hexact = loglog(data.x,data.fx,'r-','LineWidth',1.0);
for i = 1:numClusters
    plot(data.x,centroidYs(:,i));
end

set(gca,'FontSize',12);
set(gca,'Color','w');
set(gca,'XScale','log','YScale','log');
xlabel('Array Spacing (m)');
ylabel('Apparent Resistivity (\Omega-m)')

%% kmeans 2

disp('kmeans2');
numClusters = 2;
[idx,C,sumd] = kmeans(logRhoPlot',numClusters,'Distance','cityblock');
figure;
pcolor(10.^c{1},10.^c{2},N'); shading flat;
set(gca,'XScale','log','YScale','log');
hold on
for i = 1:size(C,1)
    plot(10.^C(i,:),10.^logDepthPlot(:,1),'--','LineWidth',2);
end
plot(10.^trueLogRhoPlot,10.^trueLogDepthsPlot,trueColor);
colorbar();
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'Box','on');
xlabel('Resistivity (\Omega-m)');
ylabel('Depth (m)');
legend()

figure();
subplot(3,2,1);
hist(idx);
%set(gca,'YScale','log');
title('Number of ensemble slns in each cluster');
xlabel('Cluster #');
ylabel('Frequency');

numRunsEachCluster = zeros(1,numClusters);
averageDist = zeros(1,numClusters);
centroidMisfits = zeros(1,numClusters);
centroidYs = zeros(length(data.y),numClusters);
for i = 1:numClusters
    numRunsEachCluster(i) = nnz(idx==i);
    averageDist(i) = sumd(i)/numRunsEachCluster(i);
    centroidYs(:,i) = forwardModel(xVals,10.^C(i,:),data.lambda);
    % evaluate the forward model for the maximum likelihood.;
    centroidMisfits(i) = norm(centroidYs(:,i)-data.y);
end

subplot(3,2,2);
bar(averageDist);
title('Distance from Centroid in each cluster');
ylabel('Log-distance');
xlabel('Cluster');

subplot(3,2,3);
bar(centroidMisfits)
title('Misfit')
ylabel('Data-space misfit (m)')
xlabel('Cluster');

subplot(3,2,4);
histogram(results.ensembleMisfits,100);
hold on;
yy=get(gca,'YLim');
for i = 1:numClusters
    plot(centroidMisfits(i)*[1 1],yy,'LineWidth',1);
end
set(gca,'FontSize',12);
xlabel('Log misfit');
f=gcf;
plot(medianModelMisfit*[1 1],yy,'Color',medianColor,'LineWidth',1);
plot(bestFitModelMisfit*[1 1],yy,'Color',bestFitColor,'LineWidth',1);
plot(meanModelMisfit*[1 1],yy,'Color',meanColor,'LineWidth',1);
plot(maxLikelihoodMisfit*[1 1],yy,'Color',maxLikelihoodColor,'LineWidth',1);
plot(ensembleMedianModelMisfit*[1 1],yy,ensembleMedianColor);



subplot(3,2,5);
hold on;

%Ensemble solution
for i=1:nSavedPlot
    hEnsemble=plot(xVals,yVals(:,i),'Color',ensembleColor);
end
set(gca,'Box','on');
set(gcf,'Color','w');
%True model
hdata = loglog(data.x,data.y,'r.','MarkerSize',10.0);
hexact = loglog(data.x,data.fx,'r-','LineWidth',1.0);
for i = 1:numClusters
    plot(data.x,centroidYs(:,i));
end

set(gca,'FontSize',12);
set(gca,'Color','w');
set(gca,'XScale','log','YScale','log');
xlabel('Array Spacing (m)');
ylabel('Apparent Resistivity (\Omega-m)')

%}
