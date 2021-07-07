function logRhoPlot = generateLogRhoPlot(x,depths,rhos,nzplot)

minDistL = log10(min(x));
maxDistL = log10(max(x));
zVals = logspace(minDistL,maxDistL,nzplot)'; %depth values for evaluating
numSavedRuns = size(rhos,2);

logRhoPlot = zeros(nzplot,numSavedRuns);
% will contain the (log) depths resistivities 
%of every ensemble member, formatted to all be uniform.
for i = 1:numSavedRuns %for each sln...
    logRhoPlot(:,i) = longForm(zVals,depths(:,i),rhos(:,i));
end
end