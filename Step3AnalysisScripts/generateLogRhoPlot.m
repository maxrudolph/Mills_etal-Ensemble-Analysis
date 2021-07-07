function logRhoPlot = generateLogRhoPlot(depths,rhos,zVals)
%{
Function for generating 'logRhoPlot' which is a matrix where every row is
the log-resistivity at the corresponding depth value in zVals and every
column is an ensemble member
Inputs:
    depths: depth of layer interfaces of every ensemble member (should be
        each column represents one member)
    rhos: resistivity of each layer of every ensemble member, same format
    zVals: List of depths at which to generate resistivity
%}

numSavedRuns = size(rhos,2);

logRhoPlot = zeros(length(zVals),numSavedRuns);
% will contain the (log) depths resistivities
%of every ensemble member, formatted to all be uniform.
for i = 1:numSavedRuns %for each sln...
    logRhoPlot(:,i) = longForm(zVals,depths(:,i),rhos(:,i));
end
end