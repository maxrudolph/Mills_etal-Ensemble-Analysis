function ensembleAnalysis(filename)
load(filename)

ifDebugging = true;
%% Check burn-in time and distribution of perturbation choices
if ifDebugging
    figure();
    subplot(2,1,1);
    plot(results.maxLayers);
    title('Max Layers Per Step');
    ylabel('Maximum allowed layers');
    xlabel('Step #')
    
    subplot(2,1,2);
    hist(results.allChoices);
    title('Distribution of choices');
end
    
%% Run properties: all iterations
%(likelihoods, misfits, acceptance probabilities)
    figure();
    subplot(2,2,1);
    plot(results.allLikelihoods);
    set(gca,'YScale','log');
    title('Likelihood');
    xlabel('Iteration #');
    ylabel('Likelihood');

    subplot(2,2,2);
    plot(results.allMisfits);
    title('Misfit');
    ylabel('Misfit (\Omega m)');
    xlabel('Iteration');
    
    subplot(2,2,3);
    plot(results.allProbAccepts);
    title('Acceptance Probability');
    ylabel('Probability');
    xlabel('Iteration #');
    
    
    subplot(2,2,4);
    plot(log10(results.allVars));
    title('Variance variance');
    ylabel('Variance');
    xlabel('Iteration #');
    
    %% Ensemble properties
    %(# of layers, variance, misfits)
    figure();
    subplot(1,3,1);
    hist(results.ensembleNumLayers);
    title('Number of layers');
    xlabel('Number of layers');
    ylabel('Number of saved models');

    subplot(1,3,2);
    hist(log10(results.ensembleVars));
    title('Variance');
    xlabel('log(\sigma^2)');
    ylabel('Number of saved models');
    
    subplot(1,3,3);
    plot(results.ensembleMisfits);
    title('Misfit in ensemble solutions');
    xlabel('Saved run #');
    ylabel('Misfit');
    
    %% evaluate each ensemble solution on regularly spaced grid
    nxplot=200; %resolution
    nSavedPlot = 2000; %Number of saved runs to plot
    nSaved = size(results.ensembleRhos,2);
    if nSavedPlot > nSaved %If low # of saved runs, plot all, otherwise...
        nSavedPlot = nSaved;
        ind = 1:nSavedPlot; 
    else
        ind = randperm(nSaved,nSavedPlot); 
        %... only plot some randomly chosen runs
    end
    yVals = zeros(nxplot,nSavedPlot);
    
    %%%%%%% EDIT THIS YOU FOOL %%%%%%%%%%%%%
    xVals = logspace(log10(data.x(1)),log10(data.x(end)),nxplot)';
    %%%%%%% YOU ABSOLUTE OAF %%%%%%%%%%%%%%%
    
    plotLambda = makeLambda(xVals);
    
    for i=1:nSavedPlot
        yVals(:,i) = forwardModel(squeeze(results.ensembleDepths(:,ind(i))),...
            results.ensembleRhos(:,ind(i)),plotLambda);
    end
    

    %% Data space plot: true model, data+noise, posterior solutions, etc
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
        for j = 1:size(results.ensembleDepths,2) % look at each run...
            indx = nnz(meanModelDepths(i)>=results.ensembleDepths(:,j)); 
            %find the actual layer this imaginary layer lies in
            meanModelRhos(i) = meanModelRhos(i) + ...
                log10(results.ensembleRhos(indx,j));
            % add the corresponding resistivity to the total
        end
        meanModelRhos(i) = 10.^(meanModelRhos(i)/...
            size(results.ensembleDepths,2)); %find average
    end
   
    meanModelX = logspace(minDist,maxDist,numLayersMean); 
    %new imaginary measurements
    meanModelLambda = makeLambda(meanModelX);
    
    %Maximum likelihood model
    depthPlot = zeros(length(xVals),nSaved);
rhoPlot = zeros(size(depthPlot));
for i = 1:nSaved
    %interpolate resistivity onto x-grid
    depthPlot(:,i) = log10(xVals);
    nLayer = length(~isnan(results.ensembleDepths(:,i)));
    for j = 1:nLayer
        mask = log10(xVals) >= log10(results.ensembleDepths(j,i));
        rhoPlot(mask,i) = log10(results.ensembleRhos(j,i));
    end
end
    
    % compute a bivariate histogram of resitvity values from the posterior ensemble
    [N,c]=hist3([rhoPlot(:),depthPlot(:)],{linspace(0,4,200) linspace(-1,3,100)},'CDataMode','auto');% First linspace is for log(rho), second is for log(depth)
    % at each depth, find the most likely solution (ml_rho)
    maxLikelihoodRho = zeros(size(N,2),1);
for i=1:length(xVals)
    % Use ksdensity to approximate the pdf of resistivity at this depth:
    [xi,f] = ksdensity(rhoPlot(i,:));
    [~,ind] = max(xi);
    maxLikelihoodRho(i) = 10.^f(ind);
end
    maxLikelihoodY = calculateRho1D(xVals,maxLikelihoodRho,data.lambda); % evaluate the forward model for the maximum likelihood.
    maxLikelihoodMisfit = norm(data.y-maxLikelihoodY);
    hMaxLikelihood=plot(data.x,maxLikelihoodY,'m','LineWidth',2);
    
    
    [~,ind] = sort(results.ensembleMisfits);
    medianIndex= ind(floor(length(ind)/2));
    bestIndex = ind(1);
    medianModelDepths = results.ensembleDepths(:,medianIndex);
    bestFitModelDepths = results.ensembleDepths(:,bestIndex);
    medianModelRhos = results.ensembleRhos(:,medianIndex);
    bestFitModelRhos = results.ensembleRhos(:,bestIndex);
    
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
    
legend([hexact,hdata,hMean,hMedian,hBest,hEnsembleMean,hEnsemble,...
    hMaxLikelihood],{'Exact','Data+noise','Mean model','Median model',...
    'Best fit model','Ensemble mean','Ensemble','Max Likelihood Model'},...
    'Location','southwest');

%% Model space color plot
figure();
hist3([rhoPlot(:),depthPlot(:)],[100 100],'CDataMode','auto');
[N,c] = hist3([rhoPlot(:),depthPlot(:)],[100 100],'CDataMode','auto');
pcolor(c{1},c{2},N'); shading flat;
view(2);
hold on
%trueRhoPlot = zeros(size(xVals));
%trueDepthsPlot = depthPlot(:,1);
%[trueDepths,trueRhos] = modelGen(results.maxLayers(end),thisMeasure.modelChoice);
%trueNumLayers = nnz(~isnan(trueDepths));
%for j=1:trueNumLayers
%    mask = log10(xVals) >= log10(trueDepths(j));
%    trueRhoPlot(mask) = log10(trueRhos(j));   
%end
%plot(trueRhoPlot,trueDepthsPlot,'r-','LineWidth',30);
colorbar();
set(gca,'YDir','reverse');
xlabel('log10(Resistivity)');
ylabel('log10(Depth)');
plot(log10(meanModelRhos),log10(meanModelDepths),'g-','LineWidth',2)
plot(log10
end