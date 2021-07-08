function outLRhos = longForm(zVals,inDepths,inRhos)
%{
Takes a sln with format like depths = [0,1,25,NaN,NaN]'; and rhos =
[10,390,10,NaN,NaN]'; and outputs resistivity vector of same size as zVals,
with each entry corresponding to resistivity at zVals depth. So output
would be, if
zVals = [0.1,0.102,0.104,...,0.97,1.00,1.02,...,24.9,25.4,25.9,....,1000]
then outLRhos = [10,10,10,....,10,390,390,....,390,10,10,...,10];
long form makes for more accurate plotting
%}
outLRhos = nan(size(zVals));
nLayer = nnz(~isnan(inRhos));
%...Find the number of layers in the sln...
for i = 1:nLayer %for each layer...
    mask = zVals >= inDepths(i); %Find all depths in zVals within that or lower layers
    outLRhos(mask) = log10(inRhos(i)); %set all those depths equal to next resistivity value
    %AKA Make an appropriate # of values = to that layers resistivity
end
