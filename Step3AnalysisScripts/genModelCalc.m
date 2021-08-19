function outModel = genModelCalc(inRhos,inDepths,data,color,lineStyle,...
    title,forwardModel)
%{
Generate Model From Calculation 7/7/21
    This simplifies the process of generating a 'calculated' model, ie, one
where, instead of choosing a sln from the ensemble, the resistivity values
at each depth are calculated based on some property of that depth's
resistivity among all ensemble members (like the mean or median value).
    data.y is the measured resistivity values that this model's data space
output will be compared against in order to calculate the misfit (in other
words, the 'actual data'). forwardModel is the fxn handle that will be used
to generate that data space output, and it requires data.lambda.
color, lineStyle, and title are the specs for how this model will be
plotted, they follow standard matlab conventions for the plotting
parameters Color, LineStyle, and DisplayName.
    %}
    inDepths(1) = 0; %It is assumed that the resistivity at the minimum depth
    % is the same as the resistivity at the surface.
    [shortDepths,shortRhos] = shortForm(inDepths,inRhos);
    outModel = calculatedModel(inDepths,inRhos,forwardModel(shortDepths,...
        shortRhos,data.lambda),data.y,color,lineStyle,title);
end