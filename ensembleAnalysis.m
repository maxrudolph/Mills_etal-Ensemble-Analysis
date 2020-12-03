%function ensembleAnalysis(filename)
disp('Loading data...')
load(filename,'data','forwardModel','results','measure')

ifDebugging = true;
%% 1 Figures
disp('Plotting properties...');
if ifDebugging
    things = {results.maxLayers,results.allChoices};
    histTF = [0,1,];
    titles = {'MaxLayersPerStep','Distribution of choices'};
    xlabels = {'Step #',' '};
    ylabels = {'Maximum allowed layers',' '};
    lazyPlotting(2,1,histTF,things,titles,xlabels,ylabels);   
end

%Figure 2: run properties, all iterations
things = {results.allLikelihoods, results.allMisfits,...
    results.allProbAccepts, results.allVars};
histTF = zeros(1,4);
titles = {'Likelihood','Misfit','Acceptance Probability','Variance'};
xlabels = {'Iteration #','Iteration #','Iteration #','Iteration #'};
ylabels = {'Likelihood', 'Misfit (\Omega m)','Probability','Variance'};
lazyPlotting(2,2,histTF,things,titles,xlabels,ylabels);
subplot(2,2,1)
set(gca,'YScale','log');

%Figure 3: Ensemble properties
things = {results.ensembleNumLayers, log10(results.ensembleVars),...
    results.ensembleMisfits};
histTF = [1,1,0];
titles = {'Ensemble # of layers','Ensemble Variance','Ensemble Misfit'};
xlabels = {'Number of layers','log(\sigma^2)','Saved run #'};
ylabels = {'Number of saved models','Number of saved models','Misfit'};
lazyPlotting(1,3,histTF,things,titles,xlabels,ylabels);

%% Section 2 NF: evaluate ensemble solutions on regularly spaced grid
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

%% Section 3 NF: Calculating other models
disp('Calculating models...');
%Calculated mean and calculated median
%Both calculated in log space, from the entire ensemble
nzplot=500;
zVals = 10.^linspace(log10(pBounds.depthMin),log10(pBounds.depthMax),nzplot);
numSavedRuns = size(results.ensembleRhos,2);

depthPlot = zeros(nzplot, numSavedRuns);
for i = 1:numSavedRuns
    depthPlot(:,i) = zVals;
end
%make a lot of imaginary layers
meanModelRhos = zeros(nzplot,1);
medianModelRhos = zeros(nzplot,1);

rhoPlot = zeros(size(depthPlot));
for i = 1:numSavedRuns %for each run
    nLayer = nnz(~isnan(results.ensembleDepths(:,i)));
    %...Find the number of layers in that run...
    for j = 1:nLayer
        mask = zVals >= results.ensembleDepths(j,i);
        rhoPlot(mask,i) = results.ensembleRhos(j,i);
    end
end
logDepthPlot = log10(depthPlot);
logRhoPlot = log10(rhoPlot);

for i = 1:nzplot %for each imaginary layer
    meanModelRhos(i) = 10.^(mean(logRhoPlot(i,:)));
    medianModelRhos(i) = 10.^(median(logRhoPlot(i,:)));
end

meanModelY = forwardModel(zVals,meanModelRhos,data.lambda);
medianModelY = forwardModel(zVals,medianModelRhos,data.lambda);
meanModelMisfit = norm(data.y - meanModelY);
medianModelMisfit = norm(data.y- medianModelY);

%Maximum likelihood model
%Note: I only vaguely understand what's going on here
% compute a bivariate histogram of resitvity values from the posterior ensemble
[N,c]=hist3([logRhoPlot(:),logDepthPlot(:)],...
    {linspace(-10,10,400) log10(zVals)},'CDataMode','auto');
% First linspace is for log(rho), second is for log(depth)
% at each depth, find the most likely solution (ml_rho)
maxLikelihoodRho = zeros(nzplot,1);
ks_x = linspace(log10(pBounds.rhoMin),log10(pBounds.rhoMax),1e5);
parfor i=1:nzplot
    i
    % Use ksdensity to approximate the pdf of resistivity at this depth:
    [xi,f] = ksdensity(logRhoPlot(i,:),ks_x,'bandwidth',0.01);
    [~,ind1] = max(xi);
    maxLikelihoodRho(i) = 10.^f(ind1);
end
maxLikelihoodY = forwardModel(zVals,maxLikelihoodRho,data.lambda);
% evaluate the forward model for the maximum likelihood.
maxLikelihoodMisfit = norm(data.y-maxLikelihoodY);

% Ensemble median and best fit models
[~,ind2] = sort(results.ensembleMisfits);
medianIndex= ind2(floor(length(ind2)/2));
bestIndex = ind2(1);
ensembleMedianModelDepths = results.ensembleDepths(:,medianIndex);
bestFitModelDepths = results.ensembleDepths(:,bestIndex);
ensembleMedianModelRhos = results.ensembleRhos(:,medianIndex);
bestFitModelRhos = results.ensembleRhos(:,bestIndex);
ensembleMedianModelY = forwardModel(ensembleMedianModelDepths,...
    ensembleMedianModelRhos,data.lambda);
bestFitModelY = forwardModel(bestFitModelDepths,bestFitModelRhos,...
    data.lambda);
ensembleMedianModelMisfit = results.ensembleMisfits(medianIndex);
bestFitModelMisfit = results.ensembleMisfits(bestIndex);

%% Section 4 Figure: Data space plot
disp('More plots...');
figure;
subplot(4,2,[3 5 7]);
hold on;
%pick some colors
trueColor = 'r';
medianColor = 'y';
bestFitColor = 'c';
maxLikelihoodColor = 'm';
meanColor = 'g';
ensembleMedianColor = 'k*';
ensembleColor = [200 200 200]/255;

%Ensemble solution
for i=1:nSavedPlot
    hEnsemble=plot(xVals,yVals(:,i),'Color',ensembleColor);
end
set(gca,'Box','on');
set(gcf,'Color','w');
%True model
hdata = loglog(data.x,data.y,'r.','MarkerSize',10.0);
hexact = loglog(data.x,data.fx,'r-','LineWidth',1.0);
hMean = plot(data.x,meanModelY,meanColor,'LineWidth',1);
hMedian = plot(data.x,medianModelY,medianColor,'LineWidth',1);
hBest = plot(data.x,bestFitModelY,bestFitColor,'LineWidth',1);
hMaxLikelihood=plot(data.x,maxLikelihoodY,maxLikelihoodColor,...
    'LineWidth',1);
hEnsembleMean = loglog(xVals,mean(yVals,2),'k','LineWidth',1);
hEnsembleMedian = plot(data.x,ensembleMedianModelY,ensembleMedianColor);

set(gca,'FontSize',12);
set(gca,'Color','w');
set(gca,'XScale','log','YScale','log');
xlabel('Array Spacing (m)');
ylabel('Apparent Resistivity (\Omega-m)')

legend([hexact,hdata,hEnsemble,hEnsembleMean,hEnsembleMedian,...
    hMean,hMedian,hBest,hMaxLikelihood]...
    ,{'Exact','Data+noise','Ensemble','Ensemble mean','Ensemble median',...
    'Calculated mean','Calculated median','Best fit model',...
    'Max Likelihood Model'},...
    'Position',[0.2 0.2 0.1 0.2]);
%lgd.Layout.Tile = 'south';

%% Sub-figure: Model space plot
subplot(4,2,[2 4 6 8]);
pcolor(10.^binCenters{1},10.^binCenters{2},numElements'); shading flat;
set(gca,'XScale','log','YScale','log');
%colormap(flipud(crameri('roma')));
% view(2);
hold on
trueLogRhoPlot = zeros(nzplot,1);
trueLogDepthsPlot = logDepthPlot(:,1);

[trueDepths,trueRhos] = modelGen(measure.kMax,measure.modelChoice);
trueNumLayers = nnz(~isnan(trueDepths));
for j=1:trueNumLayers
    mask = log10(zVals) >= log10(trueDepths(j));
    trueLogRhoPlot(mask) = log10(trueRhos(j));
end
plot(10.^trueLogRhoPlot,10.^trueLogDepthsPlot,trueColor);
plot(rhoPlot(:,medianIndex),10.^logDepthPlot(:,medianIndex),'Color',medianColor);
plot(rhoPlot(:,bestIndex),10.^logDepthPlot(:,bestIndex),'Color',bestFitColor);
plot(maxLikelihoodRho,zVals,'Color',maxLikelihoodColor);
plot(meanModelRhos,zVals,meanColor)
colorbar();
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'Box','on');
xlabel('Resistivity (\Omega-m)');
ylabel('Depth (m)');


%% Sub-figure: Histogram of misfit in data space
subplot(4,2,1);
histogram(results.ensembleMisfits,100);
hold on;
yy=get(gca,'YLim');
plot(medianModelMisfit*[1 1],yy,'Color',medianColor,'LineWidth',1);
plot(bestFitModelMisfit*[1 1],yy,'Color',bestFitColor,'LineWidth',1);
plot(meanModelMisfit*[1 1],yy,'Color',meanColor,'LineWidth',1);
plot(maxLikelihoodMisfit*[1 1],yy,'Color',maxLikelihoodColor,'LineWidth',1);
plot(ensembleMedianModelMisfit*[1 1],yy,ensembleMedianColor);
set(gca,'FontSize',12);
xlabel('Misfit (m)');
f=gcf;

% f.Renderer='painters';

%% Save ensemble for The Sequencer
% nsequence = 100000;
% slashpos = find(filename=='/',1,'last');
% txtfile = [filename(slashpos+1:end-3) 'csv'];
% % each row should correspond to a different resistivity profile
% csvwrite(txtfile,logRhoPlot(:,end-nsequence:end)'); % note that sequencer
% struggles with very large datasets. I choose just the last nsequence
% samples in the ensemble.

%% save the figure
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
