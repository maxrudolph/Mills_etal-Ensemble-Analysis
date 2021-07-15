function firstThreeFigures(results,saveFigures,folderName)
%{
7/14/21
Generates figures 1 to 3 based on results
%}
%Figure 1: Shows burn-in period and choice distribution (should be flat
figure('visible',~saveFigures)
subplot(2,1,1),plot(results.maxLayers)
title('MaxLayersPerStep');xlabel('Step #');ylabel('Maximum allowed layers')
subplot(2,1,2),histogram(results.allChoices,'EdgeAlpha',0);
title('Distribution of choices');
saveFigs(saveFigures,folderName,'1');

%Figure 2: run properties, all iterations
figure('visible',~saveFigures)
subjects = {results.allLikelihoods, results.allMisfits,...
    results.allProbAccepts, results.allVars};
titles = {'Likelihood','Misfit (\Omega m)','Probability','Variance'};
for i = 1:length(titles)
    subplot(2,2,i), plot(subjects{i})
    title(titles{i}),ylabel(titles{i}),xlabel('Iteration #');
end
subplot(2,2,1)
set(gca,'YScale','log');
saveFigs(saveFigures,folderName,'2');

%Figure 3: Ensemble properties
figure('visible',~saveFigures)
subplot(1,3,1),histogram(results.ensembleNumLayers,'EdgeAlpha',0);
title('# of layers'),xlabel('# of layers'),ylabel('# of saved models')
subplot(1,3,2),histogram(log10(results.ensembleVars),'EdgeAlpha',0);
title('Ensemble Variance'),xlabel('log(\sigma^2)'),ylabel('# of saved models')
subplot(1,3,3),plot(results.ensembleMisfits)
xlabel('Saved run #'),ylabel('Misfit'),title('EnsembleMisfit')
saveFigs(saveFigures,folderName,'3');
end
