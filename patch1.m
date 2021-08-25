filenames = {
    '3LayerA_0_02-Jul-2021.mat';
    '3LayerA_0.01_02-Jul-2021.mat';
    '3LayerA_0.02_02-Jul-2021.mat';
    '3LayerA_0.05_02-Jul-2021.mat';
    '3LayerA_0.1_02-Jul-2021.mat';
    '3LayerA_0.2_02-Jul-2021.mat'};

for i = 1:length(filenames)
    disp('round')
    load(['OldAnalysis_' filenames{i}]);
    load(['Ensemble_' filenames{i}],'results','data','forwardModel');
    Cdi = inv(data.Cd);
    ewre2n = zeros(1,length(results.ensembleMisfits));
    denominator = sqrt(data.y'*Cdi*data.y);
    residuals = zeros(length(data.x),length(results.ensembleMisfits));
    for j = 1:length(results.ensembleMisfits)
        residuals(:,j) = forwardModel(results.ensembleDepths(:,j),...
            results.ensembleRhos(:,j),data.lambda)-data.y;
        ewre2n(j) = sqrt(residuals(:,j)'*Cdi*residuals(:,j))/denominator;
    end
    for k = 1:length(allModels)
        allModels{k}.setWRE2N(data);
    end
    [~,ind2] = sort(ewre2n); %sort all slns by misfit
    medianIndex= ind2(floor(length(ind2)/2)); %the one with median misfit
    bestIndex = ind2(1);                      %the one with lowest misfit
    %Find their corresponding slns in the ensemble
    bestFitColor = allModels{5}.color;
    medianColor = allModels{6}.color;
    dsLineStyle = allModels{5}.lineStyle;
    dMedian = genModelInd(medianIndex,zVals,data,medianColor,dsLineStyle,...
        'DS Median',forwardModel,results);
    bestFit = genModelInd(bestIndex,zVals,data,bestFitColor,dsLineStyle,...
        'DS Best Fit',forwardModel,results);
    allModels{5} = bestFit;
    allModels{6} = dMedian;
    
        save(['Analysis_' filenames{i}],'allModels','binCenters',...
        'allClusterSets','logRhoPlot','allPartitions','numElements',...
        'nRhoBins','nzplot','xVals','yVals','zVals','residuals',...
        'ewre2n','-v7.3');
end

%end