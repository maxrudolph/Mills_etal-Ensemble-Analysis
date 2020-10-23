function ensembleAnalysis(filename)
load(filename)

ifDebugging = true;

if ifDebugging
    figure();
    subplot(2,1,1);
    plot(results.maxLayers);
    title('Max Layers Per Step');
    ylabel('Maximum allowed layers');
    xlabel('Step #')
    
    subplot(2,1,2);
    hist(results.storedChoices);
    title('Distribution of choices');
end
    
%% Run properties
    figure();
    %figure 1: all iterations, properties changing over time
    subplot(1,3,1);
    plot(results.storedLikelihoods);
    set(gca,'YScale','log');
    title('Likelihood');
    xlabel('Iteration #');
    ylabel('Likelihood');

    subplot(1,3,2);
    plot(results.allMisfits);
    title('Misfit');
    ylabel('Misfit (\Omega m)');
    xlabel('Iteration');
    
    subplot(1,3,3);
    plot(results.storedProbAccepts);
    title('Acceptance Probability');
    ylabel('Probability');
    xlabel('Iteration #');
    
    %% Plot the ensemble properties vs. iteration
    figure();
    subplot(1,3,1);
    hist(results.numLayers);
    title('Number of layers');
    xlabel('Number of layers');
    ylabel('Number of saved models');

    subplot(1,3,2);
    hist(log10(results.storedSavedVars));
    title('Variance');
    xlabel('log(\sigma^2)');
    ylabel('Number of saved models');
    
    subplot(1,3,3);
    plot(results.ensembleMisfits);
    title('Misfit in ensemble solutions');
    xlabel('Saved run #');
    ylabel('Misfit');
    
    %% evaluate the solution for each ensemble member on regularly spaced grid
    nxplot=200; %resolution
    nSavedPlot = 2000; %Number of saved runs to plot
    nSaved = size(results.storedRhos,2);
    if nSavedPlot > nSaved %If low # of saved runs, plot all, otherwise...
        nSavedPlot = nSaved;
        ind = 1:nSavedPlot; 
    else
        ind = randperm(nSaved,nSavedPlot); 
        %... only plot some randomly chosen runs
    end
    yVals = zeros(nxplot,nSavedPlot);
    xVals = logspace(log10(data.x(1)),log10(data.x(end)),nxplot)';
    %Needs further modularization.... :-(
    plotLambda = makeLambda(xVals);
    
    for i=1:nSavedPlot
        yVals(:,i) = forwardModel(squeeze(results.storedDepths(:,ind(i))),...
            results.storedRhos(:,ind(i)),plotLambda);
    end
    

    %% Plot true model, the data+noise, the posterior solutions, and mean model
    figure;
    subplot(4,2,[3 5 7]);
    hold on;
    
    %Ensemble solution    
    for i=1:nSavedPlot
        hEnsemble=plot(xVals,yVals(:,i),'Color',[200 200 200]/255); 
    end
    set(gca,'Box','on');
    set(gcf,'Color','w');
    %True model
    hdata = loglog(data.x,data.y,'r.','MarkerSize',10.0);
    hexact = loglog(data.x,data.fx,'r-','LineWidth',1.5);
    
    %Mean model
    minDist = log10(data.x(1));
    maxDist = log10(data.x(end));
    
    numLayersMean = results.maxLayers(end)*10;
    meanModelDepths = [0,logspace(0,maxDist, numLayersMean-1)];
    meanModelRhos = zeros(numLayersMean,1);
    for i = 1:size(meanModelDepths,2) %for each imaginary layer...
        for j = 1:size(results.storedDepths,2) % look at each run...
            indx = nnz(meanModelDepths(i)>=results.storedDepths(:,j)); 
            %find the actual layer this imaginary layer lies in
            meanModelRhos(i) = meanModelRhos(i) + ...
                log10(results.storedRhos(indx,j));
            % add the corresponding resistivity to the total
        end
        meanModelRhos(i) = 10.^(meanModelRhos(i)/...
            size(results.storedDepths,2)); %find average
    end
   
    meanModelX = logspace(minDist,maxDist,numLayersMean); 
    %new imaginary measurements
    meanModelLambda = makeLambda(meanModelX);
    
    [~,ind] = sort(results.ensembleMisfits);
    medianIndex= ind(floor(length(ind)/2));
    bestIndex = ind(1);
    medianModelDepths = results.storedDepths(:,medianIndex);
    bestFitModelDepths = results.storedDepths(:,bestIndex);
    medianModelRhos = results.storedRhos(:,medianIndex);
    bestFitModelRhos = results.storedRhos(:,bestIndex);
    
    hMean = plot(meanModelX,forwardModel(meanModelDepths,meanModelRhos,...
        meanModelLambda),'g-','LineWidth',2);
    hMedian = plot(data.x,forwardModel(medianModelDepths,medianModelRhos,...
        data.lambda),'LineWidth',2);
    hBest = plot(data.x,forwardModel(bestFitModelDepths,...
        bestFitModelRhos,data.lambda),'LineWidth',2);
    
    hEnsembleMean = loglog(xVals,mean(yVals,2),'k','LineWidth',1);
    
    set(gca,'FontSize',12);
    set(gca,'Color','w');
    set(gca,'XScale','log');
    set(gca,'YScale','log');

    xlabel('Array spacing (m)');
    ylabel('Apparent Resistivity (\Omega-m)')
    
legend([hexact,hdata,hMean,hMedian,hBest,hEnsembleMean,hEnsemble],...
    {'Exact','Data+noise','Mean model','Median model','Best fit model',...
    'Ensemble mean','Ensemble'},'Location','southwest');

%% Color plot
depthPlot = zeros(length(xVals),nSaved);
rhoPlot = zeros(size(depthPlot));
for i = 1:nSaved
    %interpolate resistivity onto x-grid
    depthPlot(:,i) = log10(xVals);
    nLayer = length(~isnan(results.storedDepths(:,i)));
    for j = 1:nLayer
        mask = log10(xVals) >= log10(results.storedDepths(j,i));
        rhoPlot(mask,i) = log10(results.storedRhos(j,i));
    end
end

figure();
hist3([rhoPlot(:),depthPlot(:)],[100 100],'CDataMode','auto');
view(2);
hold on
trueRhoPlot = zeros(size(xVals));
trueDepthsPlot = depthPlot(:,1);
[trueDepths,trueRhos] = modelGen(results.maxLayers(end),thisMeasure.modelChoice);
trueNumLayers = nnz(~isnan(trueDepths));
for j=1:trueNumLayers
    mask = log10(xVals) >= log10(trueDepths(j));
    trueRhoPlot(mask) = log10(trueRhos(j));   
end
plot(trueRhoPlot,trueDepthsPlot,'r-','LineWidth',30);
colorbar();
set(gca,'YDir','reverse');
xlabel('log10(Resistivity)');
ylabel('log10(Depth)');
plot(log10(meanModelRhos),log10(meanModelDepths),'g-','LineWidth',30)
end