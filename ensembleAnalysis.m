function ensembleAnalysis(results,data,model)
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
        yVals(:,i) = model(squeeze(results.storedDepths(:,ind(i))),...
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
    
    %{
    % find solution with median misfit from ensemble
    [meanModel, medianModel, bestFitModel] = createMeanModel(...
        length(results.storedRhos(:,1)),measure,results);
    
   
    hMedian = plot(data.x,model(medianModel.depths,medianModel.rhos,...
        data.lambda),'LineWidth',2);
    hBest = plot(data.x,model(bestFitModel.depths,bestFitModel.rhos,...
        data.lambda),'LineWidth',2);

hmean=loglog(x_vals,mean(y_vals,2),'k','LineWidth',1); %Ensemble mean


hEnsembleMean = loglog(newX,calculateRho1D(meanDepths,meanRhos,meanLambda),'g','LineWidth',2);
%calculated mean
set(gca,'FontSize',12);
set(gca,'Color','w');
set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('Array spacing (m)');
ylabel('Apparent Resistivity (\Omega-m)')
    
hmean=loglog(xVals,mean(yVals,2),'k','LineWidth',1); %Ensemble mean
hOtherMean = loglog(newX,calculateRho1D(meanDepths,meanRhos,meanLambda),'g*');
%calculated mean
set(gca,'FontSize',16);
set(gca,'Color','w');
set(gca,'XScale','log');
set(gca,'YScale','log');
legend([hexact,hdata,h1,hmean,hOtherMean],...
    {'Exact','Data+noise','Ensemble','Ensemble mean','Calculated mean'},...
    'Location','southwest');

%% Plot the model
depth_plot = zeros(length(xVals),nSaved);
rho_plot = zeros(size(depth_plot));
for i=1:nSaved
    % interpolate resistivity onto x-grid
    depth_plot(:,i) = log10(xVals);
    nlayer = length(~isnan(results.storedDepths(:,i)));
    for j=1:nlayer
        mask = log10(xVals) >= log10(results.storedDepths(j,i));
        rho_plot(mask,i) = log10(results.storedRhos(j,i));
    end
end
figure();
hist3([rho_plot(:),depth_plot(:)],[100 100],'CDataMode','auto');
view(2);
hold on
true_rho_plot = zeros(size(xVals));
true_depths_plot = depth_plot(:,1);
for j=1:trueNumLayers
    mask = log10(xVals) >= log10(trueDepths(j));
    true_rho_plot(mask) = log10(trueRhos(j));   
end
plot(true_rho_plot,true_depths_plot,'r');
colorbar();
set(gca,'YDir','reverse');
xlabel('log10(Resistivity)');
ylabel('log10(Depth)');
plot(log10(meanRhos),log10(meanDepths),'g')
    %}
end