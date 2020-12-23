%function ensembleAnalysis(filename)
rng(1);
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

dMean = makeCalculatedModel(zVals,10.^(mean(logRhoPlot,2)),data,...
    forwardModel,'g','Data Space Mean');
dMedian = makeCalculatedModel(zVals,10.^(median(logRhoPlot,2)),data,...
    forwardModel,'y','Data Space Median');

% Ensemble median and best fit models
[~,ind2] = sort(results.ensembleMisfits);
medianIndex= ind2(floor(length(ind2)/2));
bestIndex = ind2(1);
mMedian = makeCalculatedModel(results.ensembleDepths(:,medianIndex),...
    results.ensembleRhos(:,medianIndex),data,forwardModel,'k',...
    'Model Space Median');
bestFit = makeCalculatedModel(results.ensembleDepths(:,bestIndex),...
    results.ensembleRhos(:,bestIndex),data,forwardModel,'c',...
    'Best Fit Model');
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
    % Use ksdensity to approximate the pdf of resistivity at this depth:
    [pdfYVals,pdfXVals] = ksdensity(logRhoPlot(:,i),ksRho,'bandwidth',.05);
    [~,ind1] = max(pdfYVals);
    maxLikelihoodRho(i) = 10.^pdfXVals(ind1);
end
logRhoPlot = logRhoPlot';
maxLikelihood = makeCalculatedModel(zVals,maxLikelihoodRho,data,...
    forwardModel,'m','Maximum Likelihood Model');

%Setup true model
[trueDepths,trueRhos] = modelGen(measure.kMax,measure.modelChoice);
trueModel = makeCalculatedModel(trueDepths,trueRhos,data,forwardModel,...
    'r','True Model');

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
allModels = {trueModel,dMean,dMedian,bestFit,maxLikelihood,...
    mMedian};
plots = cell(1,6);
for iPlot = 1:6
    plots{iPlot} = loglog(data.x,allModels{iPlot}.y,'-','Color',...
        allModels{iPlot}.color);
end

hdata = loglog(data.x,data.y,'r.','MarkerSize',10.0);
hmMean = loglog(xVals,mean(yVals,2),'ks');

legend([plots{1},plots{2},plots{3},plots{4},plots{5},plots{6},hdata,...
    hmMean,hEnsemble],{allModels{1}.displayName,allModels{2}.displayName,...
    allModels{3}.displayName,allModels{4}.displayName,...
    allModels{5}.displayName,allModels{6}.displayName,'Data+noise',...
    'Model Space mean','Ensemble'},'Position',[0.2 0.2 0.1 0.2]);

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
    'Color',mMedian.color);
plot(rhoPlot(:,bestIndex),10.^logDepthPlot(:,bestIndex),'-',...
    'Color',bestFit.color);
plot(maxLikelihood.rhos,maxLikelihood.depths,'-',...
    'Color',maxLikelihood.color);
plot(dMean.rhos,dMean.depths,'-','Color',dMean.color);
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
for iPlot = 2:6
    plot(allModels{iPlot}.misfit*[1 1],yy,'-','Color',allModels{iPlot}.color);
end
set(gca,'FontSize',12);
xlabel('Misfit (m)');

%% GM Model
disp('Next phase...');
rng(1);
%downsample
downsampleNumber = floor(size(results.ensembleRhos,2))%/10);
gmPlotIndex = randperm(size(results.ensembleRhos,2),downsampleNumber);

%Step one: use k-means clustering
numClusters = 6;
clust = zeros(downsampleNumber,numClusters);
gmRhoPlot = logRhoPlot(:,gmPlotIndex);
parfor i=1:numClusters
    disp(i);
    clust(:,i) = kmeans(gmRhoPlot',i);
end

eva = evalclusters(gmRhoPlot',clust,'CalinskiHarabasz')
options = statset('Display','final');
GMModel = fitgmdist(gmRhoPlot',eva.OptimalK,...
    'Options',options,'RegularizationValue',1e-6);
%AIC = zeros(1,numClusters);
%GMModels = cell(1,numClusters);
%options = statset('MaxIter',500);
%parfor k = 1:numClusters
%    disp(k)
%    GMModels{k} = fitgmdist(logRhoPlot',k,'Options',options,'CovarianceType','diagonal');
%    AIC(k)= GMModels{k}.AIC;
%end

%[minAIC,numComponents] = min(AIC);
%numComponents;
%BestModel = GMModels{numComponents};

figure()
pcolor(10.^binCenters{1},10.^binCenters{2},numElements'); shading flat;
set(gca,'XScale','log','YScale','log');
hold on
plot(10.^GMModel.mu,10.^logDepthPlot(:,1)','--','DisplayName','GM means');
plot(10.^trueLogRhoPlot,10.^trueLogDepthsPlot,trueModel.color,'DisplayName','True model');
colorbar();
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'Box','on');
xlabel('Resistivity (\Omega-m)');
ylabel('Depth (m)');
legend()
title('GM Stuff')
figure()
histogram(results.ensembleMisfits,100);
hold on;
yy=get(gca,'YLim');
gmYs = zeros(size(data.y,1),eva.OptimalK);
gmMisfits = zeros(1,eva.OptimalK);
for i = 1:eva.OptimalK
    gmYs(:,i) = forwardModel(zVals,10.^GMModel.mu(i,:)',data.lambda);
    % evaluate the forward model for the maximum likelihood.;
    gmMisfits(i) = norm(gmYs(:,i)-data.y);
    plot(gmMisfits(i)*[1 1],yy,'--');
end

for iPlot = 2:6
    plot(allModels{iPlot}.misfit*[1 1],yy,'-','Color',allModels{iPlot}.color);
end

%%
figure()
[idx,C,sumd] = kmeans(logRhoPlot',eva.OptimalK);
pcolor(10.^binCenters{1},10.^binCenters{2},numElements'); shading flat;
set(gca,'XScale','log','YScale','log');
hold on
for i = 1:size(C,1)
    plot(10.^C(i,:),10.^logDepthPlot(:,1),'--','LineWidth',2,'DisplayName',...
        'Centroids');
end
plot(10.^trueLogRhoPlot,10.^trueLogDepthsPlot,trueModel.color,'DisplayName',...
    'True model');
colorbar();
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'Box','on');
xlabel('Resistivity (\Omega-m)');
ylabel('Depth (m)');
legend()

figure();
subplot(3,2,1);
histogram(idx);
%set(gca,'YScale','log');
title('Number of ensemble slns');
xlabel('Cluster #');
ylabel('Frequency');

numRunsEachCluster = zeros(1,eva.OptimalK);
averageDist = zeros(1,eva.OptimalK);
centroidMisfits = zeros(1,eva.OptimalK);
centroidYs = zeros(length(data.y),eva.OptimalK);
for i = 1:eva.OptimalK
    numRunsEachCluster(i) = nnz(idx==i);
    averageDist(i) = sumd(i)/numRunsEachCluster(i);
    centroidYs(:,i) = forwardModel(zVals,10.^C(i,:),data.lambda);
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
for i = 1:eva.OptimalK
    plot(centroidMisfits(i)*[1 1],yy,'--','LineWidth',1,'DisplayName',...
        'Centroids');
end
set(gca,'FontSize',12);
xlabel('Log misfit');
f=gcf;
for iPlot = 2:6
    plot(allModels{iPlot}.misfit*[1 1],yy,'-','Color',allModels{iPlot}.color,...
        'DisplayName',allModels{iPlot}.displayName);
end

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
for i = 1:eva.OptimalK
    plot(data.x,centroidYs(:,i));
end

set(gca,'FontSize',12);
set(gca,'Color','w');
set(gca,'XScale','log','YScale','log');
xlabel('Array Spacing (m)');
ylabel('Apparent Resistivity (\Omega-m)')

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


%% AIC and BIC vs regularization value
disp('testing phase...');
rng(1);
%downsample
downsampleNumber = floor(size(results.ensembleRhos,2)/10);
gmPlotIndex = randperm(size(results.ensembleRhos,2),downsampleNumber);

%Step one: use k-means clustering
numClusters = 6;
clust = zeros(downsampleNumber,numClusters);
gmRhoPlot = logRhoPlot(:,gmPlotIndex);
parfor i=1:numClusters
    disp(i);
    clust(:,i) = kmeans(gmRhoPlot',i);
end

AIC = zeros(1,10);
BIC = zeros(1,10);

for iReg = 1:10
    eva = evalclusters(gmRhoPlot',clust,'CalinskiHarabasz');
    options = statset('Display','final');
    GMModel = fitgmdist(gmRhoPlot',eva.OptimalK,...
        'Options',options,'RegularizationValue',1*(10^(-iReg)));
    AIC(iReg) = GMModel.AIC;
    BIC(iReg) = GMModel.BIC;
    disp(iReg)
end

figure;
plot([1:10],AIC,'*','DisplayName','AIC');
hold on
plot([1:10],BIC,'^','DisplayName','BIC');
xlabel('Regularization Value = 1e(-x)')
ylabel('AIC or BIC');
legend;
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
%}
%% k-means
%{
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
histogram(idx);
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
plot(mMedianModelMisfit*[1 1],yy,mMedianColor);


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
plot(mMedianModelMisfit*[1 1],yy,mMedianColor);



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

function plotTrueModelSpace()
end