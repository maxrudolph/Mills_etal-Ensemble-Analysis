function ensembleAnalysis(results,measure,model)
    %% Run properties
    figure();
    subplot(2,1,1);
    plot(results.storedLikelihoods);
    set(gca,'YScale','log');
    title('Likelihood');
    xlabel('Iteration #');
    ylabel('Likelihood');

    subplot(2,1,2);
    plot(results.allMisfits);
    title('Misfit');
    ylabel('Misfit (\Omega m)');
    xlabel('Iteration');
    
    %% Plot the ensemble properties vs. iteration
    figure();
    subplot(2,1,1);
    hist(results.numLayers);
    title('Number of layers');
    xlabel('Number of layers');
    ylabel('Number of saved models');

    subplot(2,1,2);
    hist(log10(results.storedSavedVars));
    title('Variance');
    xlabel('log(\sigma^2)');
    ylabel('Number of saved models');
    
    
    %% evaluate the solution for each ensemble member on regularly spaced grid
    nxplot=200;
    nSavedPlot = 2000;
    nSaved = size(results.storedRhos,2);
    if nSavedPlot > nSaved
        nSavedPlot = nSaved;
        ind = 1:nSavedPlot;
    else
        ind = randperm(nSaved,nSavedPlot);
    end
    yVals = zeros(nxplot,nSavedPlot);
    xVals = logspace(measure.minDist,measure.maxDist,nxplot)';
    plotLambda = makeLambda(xVals);
    
    for i=1:nSavedPlot
        yVals(:,i) = model(squeeze(results.storedDepths(:,ind(i))),...
            results.storedRhos(:,ind(i)),plotLambda);
    end
    

    %% Plot true model, the data+noise, the posterior solutions, and mean model
    figure;
    subplot(4,2,[3 5 7]);
    hold on;
    for i=1:nSavedPlot
        hEnsemble=plot(xVals,yVals(:,i),'Color',[200 200 200]/255); %The ensemble solutions
    end
    % find solution with median misfit from ensemble
    [meanModel, medianModel, bestFitModel] = createMeanModel(...
        length(results.storedRhos(:,1)),measure,results);
    
   
    hMedian = plot(measure.x,model(medianModel.depths,medianModel.rhos,...
        measure.lambda),'LineWidth',2);
    hBest = plot(measure.x,model(bestFitModel.depths,bestFitModel.rhos,...
        measure.lambda),'LineWidth',2);

   %{
set(gca,'Box','on');
set(gcf,'Color','w');
hdata=loglog(x,y,'r.','MarkerSize',10.0); %The measurements
hmean=loglog(x_vals,mean(y_vals,2),'k','LineWidth',1); %Ensemble mean
hexact=loglog(x,fx,'r-','LineWidth',1.5); %The true model

hEnsembleMean = loglog(newX,calculateRho1D(meanDepths,meanRhos,meanLambda),'g','LineWidth',2);
%calculated mean
set(gca,'FontSize',12);
set(gca,'Color','w');
set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('Array spacing (m)');
ylabel('Apparent Resistivity (\Omega-m)')
    
    
    
    
    
    
    
    %% Plot true model, the data+noise, the posterior solutions, and mean model
    figure;
    hold on;
    for i=1:nSaved
        h1=plot(xVals,yVals(:,i),'Color',[200 200 200]/255); %The ensemble solutions
    end
    set(gca,'Box','on');
    set(gcf,'Color','w');
    hdata=loglog(x,y,'r.','MarkerSize',10.0); %The measurements
hmean=loglog(xVals,mean(yVals,2),'k','LineWidth',1); %Ensemble mean
hexact=loglog(x,fx,'r-','LineWidth',1.5); %The true model
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