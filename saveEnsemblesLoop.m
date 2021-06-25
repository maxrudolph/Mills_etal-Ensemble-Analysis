function saveEnsemblesLoop
noiseLevels = [0,0.01,0.02,0.05,0.1,0.2];
noiseLevels = repmat(noiseLevels,1,3)';
subStructs = {'1LayerA','3LayerA','4LayerA'};
subStructs = repmat(subStructs,6,1);
subStructs = vertcat(subStructs{:,1},subStructs{:,2},subStructs{:,3});



%Create data files
tic
parfor i = 1:length(noiseLevels)+3
    if i>length(noiseLevels)
        a = createSyntheticData(0,'subStructChoice',...
            subStructs((i-length(noiseLevels))*length(subStructs)/3,:))
        b = inversion(a,'priorOn',true);
    else %do the priors for each subStruct
        a = createSyntheticData(noiseLevels(i),'subStructChoice',subStructs(i,:));
        b = inversion(a);
    end
end
toc

