function outModel = genModelInd(index,inDepths,data,color,lineStyle,title,...
    forwardModel,results,options,pBounds)
%{
Generate Model From Index 7/7/2021
    Some models, like one with best fit to the data in data space, are not
calculated from all ensemble members, but are just an actual member chosen
from the ensemble. Pass the results structure (the entire library of
ensemble member properties) and the index of the desired sln, and it will
assemble the calculatedModel, including converting it to long form for
plotting. color, lineStyle, and title follow standard matlab conventions
for Color, LineStyle, and DisplayName when plotting. forwardModel is the
fxn handle used to find the data space output of inDepths, it requires
data.lambda in order to be evaluated at the same points as the 'actual
data' were. data.y (the 'actual data') are also required for comparing
values and calculating misfit. And inDepths is the long array of
(preferably log-spaced) to generate the resistivities for plotting.
%}
rhos = results.ensembleRhos(:,index);
depths = results.ensembleDepths(:,index);
if options.piecewiseLinear
    [ztmp,rhotmp] = piecewiseLinearSolution(depths,rhos,pBounds,inDepths);
    rhoPlot = rhotmp;
else
    rhoPlot = longForm(inDepths,depths,rhos);
end


outModel = calculatedModel(inDepths,rhoPlot,...
    forwardModel(depths,rhos,data.lambda),data.y,color,lineStyle,title);
outModel.setWRE2N(data);
end