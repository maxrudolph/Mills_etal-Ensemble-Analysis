%function ensembleAnalysis(filename)
load(filename)


ifDebugging = true;
%% Figure: Check burn-in time and distribution of perturbation choices
if ifDebugging
    figure();
    subplot(2,1,1);
    plot(results.maxLayers);
    title('Max Layers Per Step');
    ylabel('Maximum allowed layers');
    xlabel('Step #')
    
    subplot(2,1,2);
    histogram(results.allChoices);
    title('Distribution of choices');
end
    
%% Figure: Run properties, all iterations
%(likelihoods, misfits, acceptance probabilities)
disp('Plotting properties...');
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
    
    %% Figure: Ensemble properties
    %(# of layers, variance, misfits)
    figure();
    subplot(1,3,1);
    histogram(results.ensembleNumLayers);
    title('Number of layers');
    xlabel('Number of layers');
    ylabel('Number of saved models');

    subplot(1,3,2);
    histogram(log10(results.ensembleVars));
    title('Variance');
    xlabel('log(\sigma^2)');
    ylabel('Number of saved models');
    
    subplot(1,3,3);
    %f=fit(1:size(results.ensembleMisfits,2),results.ensembleMisfits,'poly1')
    plot(results.ensembleMisfits);
    %plot(f);
    title('Misfit in ensemble solutions');
    xlabel('Saved run #');
    ylabel('Misfit');
    
    
    %% NF: evaluate each ensemble solution on regularly spaced grid
    disp('Evaluating ensemble...');
    minDist = log10(data.x(1)); %needed later
    maxDist = log10(data.x(end));    
    nxplot=200; %resolution of measurement points
    nSavedPlot = 2000; %Number of saved runs to plot
    if nSavedPlot > size(results.ensembleRhos,2) 
        %If low # of saved runs, plot all, otherwise...
        nSavedPlot = size(results.ensembleRhos,2);
        runPlotIndex = 1:nSavedPlot; 
    else %... only plot a random subset of them
        runPlotIndex = randperm(size(results.ensembleRhos,2),nSavedPlot); 
    end
    yVals = zeros(nxplot,nSavedPlot);
    xVals = logspace(minDist,maxDist,nxplot)';
    lambdaForXVals = makeLambda(xVals);
    
    for i=1:nSavedPlot
        yVals(:,i) = forwardModel(squeeze(results.ensembleDepths(:,runPlotIndex(i))),...
            results.ensembleRhos(:,runPlotIndex(i)),lambdaForXVals);
    end
    
    %% NF: Calculating other models
    disp('Calculating models...');
    %Mean model and calculated median model
    %Both calculated in log space, from the entire ensemble
    numSavedRuns = size(results.ensembleRhos,2);
    depthPlot = zeros(length(xVals), numSavedRuns);
    for i = 1:numSavedRuns
        depthPlot(:,i) = xVals;
    end
    %make a lot of imaginary layers
    meanModelRhos = zeros(nxplot,1);
    medianModelRhos = zeros(nxplot,1);
    
    rhoPlot = zeros(size(depthPlot));
    for i = 1:numSavedRuns %for each run
        nLayer = nnz(~isnan(results.ensembleDepths(:,i)));
        %Find the number of layers in that run
        for j = 1:nLayer
            mask = xVals >= results.ensembleDepths(j,i);
            rhoPlot(mask,i) = results.ensembleRhos(j,i);
        end
    end
    logDepthPlot = log10(depthPlot);
    logRhoPlot = log10(rhoPlot);
    
    for i = 1:nxplot %for each imaginary layer
        meanModelRhos(i) = 10.^(mean(logRhoPlot(i,:)));
        medianModelRhos(i) = 10.^(median(logRhoPlot(i,:)));
    end

    meanModelY = forwardModel(xVals,meanModelRhos,data.lambda);
    medianModelY = forwardModel(xVals,medianModelRhos,data.lambda);
    meanModelMisfit = norm(data.y - meanModelY);
    medianModelMisfit = norm(data.y- medianModelY);
    assert(all(size(data.y)==size(meanModelY)));
    assert(all(size(data.y)==size(medianModelY)))
    
    %Maximum likelihood model    
    %Note: I only vaguely understand what's going on here
    % compute a bivariate histogram of resitvity values from the posterior ensemble
    [N,c]=hist3([logRhoPlot(:),logDepthPlot(:)],...
        {linspace(-10,10,400) linspace(minDist,maxDist,nxplot)},'CDataMode','auto');
    % First linspace is for log(rho), second is for log(depth)
    % at each depth, find the most likely solution (ml_rho)
    maxLikelihoodRho = zeros(size(N,2),1);
    for i=1:length(xVals)
        % Use ksdensity to approximate the pdf of resistivity at this depth:
        [xi,f] = ksdensity(logRhoPlot(i,:));
        [~,ind1] = max(xi);
        maxLikelihoodRho(i) = 10.^f(ind1);
    end
    maxLikelihoodY = forwardModel(xVals,maxLikelihoodRho,data.lambda); 
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
    ensembleMedianModelMisfit = norm(data.y - ensembleMedianModelY);
    bestFitModelMisfit = norm(data.y - bestFitModelY);
    
    %% Figure: Data space plot
    disp('Last plot...');
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
    pcolor(10.^c{1},10.^c{2},N'); shading flat;
    set(gca,'XScale','log','YScale','log');
    %colormap(flipud(crameri('roma')));
    % view(2);
    hold on
    trueRhoPlot = zeros(size(xVals));
    trueDepthsPlot = depthPlot(:,1);
    [trueDepths,trueRhos] = modelGen(options.kMax,measure.modelChoice);
    trueNumLayers = nnz(~isnan(trueDepths));
    for j=1:trueNumLayers
        mask = log10(xVals) >= log10(trueDepths(j));
        trueRhoPlot(mask) = log10(trueRhos(j));
    end
    plot(10.^trueRhoPlot,10.^trueDepthsPlot,trueColor);
    plot(10.^rhoPlot(:,medianIndex),10.^depthPlot(:,medianIndex),'Color',medianColor);
    plot(10.^rhoPlot(:,bestIndex),10.^depthPlot(:,bestIndex),'Color',bestFitColor);
    plot(maxLikelihoodRho,xVals,'Color',maxLikelihoodColor);   
    plot(meanModelRhos,xVals,meanColor)
    colorbar();
    set(gca,'YDir','reverse');
    set(gca,'FontSize',12);
    set(gca,'Box','on');
    xlabel('Resistivity (\Omega-m)');
    ylabel('Depth (m)');

%% Sub-figure: Histogram of misfit
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
    f.Renderer='painters';
%% save the figure
% export_fig([results_file(1:end-4) '_ensemble.eps'])
%exportgraphics(gcf,[results_file(1:end-4) '_ensemble.eps'],'ContentType','vector');

%end