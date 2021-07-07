function outLRhos = longForm(zVals,inDepths,inRhos)

outLRhos = nan(size(zVals));
nLayer = nnz(~isnan(inRhos));
    %...Find the number of layers in that sln...
    for i = 1:nLayer %for each layer...
        mask = zVals >= inDepths(i);
        outLRhos(mask) = log10(inRhos(i));
        %Make an appropriate # of values = to that layers resistivity
    end
    