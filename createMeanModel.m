function [meanModel,medianModel,bestFitModel] = ...
    createMeanModel(numLayersMean,measure,results)
    %Model represents the mean of all models (at each depth, the mean
    %resistivity across all saved models for that depth)
    meanModel.depths = [0,logspace(0,measure.maxDist,numLayersMean-1)];
    meanModel.rhos = zeros(numLayersMean,1);
    for i = 1:size(meanModel.depths,2) %for each imaginary layer...
        for j = 1:size(results.storedDepths,2) % look at each run...
            indx = nnz(meanModel.depths(i)>=results.storedDepths(:,j)); 
            %find the actual layer this imaginary layer lies in
            meanModel.rhos(i) = meanModel.rhos(i) + ...
                log10(results.storedRhos(indx,j));
            % add the corresponding resistivity to the total
        end
        meanModel.rhos(i) = 10.^(meanModel.rhos(i)/...
            size(results.storedDepths,2)); %find average
    end
    
    newX = logspace(measure.minDist,measure.maxDist,numLayersMean); 
    %new imaginary measurements
    meanModel.lambda = makeLambda(newX);
    
    %% Median, best
    % find solution with median misfit from ensemble
    [~,ind] = sort(results.ensembleMisfits);
    medianIndex= ind(floor(length(ind)/2));
    bestIndex = ind(1);
    medianModel.depths = results.storedDepths(:,medianIndex);
    bestFitModel.depths = results.storedDepths(:,bestIndex);
    medianModel.rhos = results.storedRhos(:,medianIndex);
    bestFitModel.rhos = results.storedRhos(:,bestIndex);
    
    %medianModel = plot(x,model(results.storedDepths(:,medianIndex),...
     %   results.storedRhos(:,medianIndex),measure.lambda),'LineWidth',2);
    %bestFitModel = plot(x,model(results.storedDepths(:,bestIndex),...
     %   results.storedRhos(:,bestIndex),measure.lambda),'LineWidth',2);
end