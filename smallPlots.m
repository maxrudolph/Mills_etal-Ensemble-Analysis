function smallPlots(results,saveFigures,folderName,visibility)
%Figure 1: debugging stuff
things = {results.maxLayers,results.allChoices};
histTF = [0,1,];
titles = {'MaxLayersPerStep','Distribution of choices'};
xlabels = {'Step #',' '};
ylabels = {'Maximum allowed layers',' '};
lazyPlotting(2,1,histTF,things,titles,xlabels,ylabels,visibility);
saveFigs(saveFigures,folderName,'1');

%Figure 2: run properties, all iterations
things = {results.allLikelihoods, results.allMisfits,...
    results.allProbAccepts, results.allVars};
histTF = zeros(1,4);
titles = {'Likelihood','Misfit','Acceptance Probability','Variance'};
xlabels = {'Iteration #','Iteration #','Iteration #','Iteration #'};
ylabels = {'Likelihood', 'Misfit (\Omega m)','Probability','Variance'};
lazyPlotting(2,2,histTF,things,titles,xlabels,ylabels,visibility);
subplot(2,2,1)
set(gca,'YScale','log');
saveFigs(saveFigures,folderName,'2');


%Figure 3: Ensemble properties
things = {results.ensembleNumLayers, log10(results.ensembleVars),...
    results.ensembleMisfits};
histTF = [1,1,0];
titles = {'Ensemble # of layers','Ensemble Variance','Ensemble Misfit'};
xlabels = {'Number of layers','log(\sigma^2)','Saved run #'};
ylabels = {'Number of saved models','Number of saved models','Misfit'};
lazyPlotting(1,3,histTF,things,titles,xlabels,ylabels,visibility);
saveFigs(saveFigures,folderName,'3');
end
