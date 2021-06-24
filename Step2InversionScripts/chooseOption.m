function output = chooseOption(currentLayers,maxLayers,numOptions)
%{
6/21/2021
A subscript for mcmcAlgorithm. This script is called at each step in the
main loop. It chooses the way the accepted sln will be edited to form the
proposed sln: either adding or removing a layer, changing the depth of an
interface, changing the resistivity of a layer, or changing the variance
(if that is an allowed option). The probabilities and availabilities of
different options need to be controlled depending on the state of the
current accepted Sln (for instance, if the accepted sln is only one layer,
then deleteLayer should not be available, but also, the probability to add
a new layer should not be larger than it normally is, or that would bias
the program toward multi-layer models).
%}
if currentLayers == 1 %If one layer model...
    % prob(deleteLayer) = prob(perturbDepth) = 0;
    % If numOptions = 4, prob(addLayer) = 25%, prob(perturbRho) = 75%;
    % if numOptions = 5, prob(addLayer) = 20%, prob(perturbRho) = 40%,
    % prob(changeVar) = 40%;
    tmp = rand(); %uniformly distributed
    if tmp<(1/numOptions)
        output = 3; %Add new layer.
    else
        availableOptions = [4:numOptions];
        output = availableOptions(randi(length(availableOptions)));
    end
elseif currentLayers >= maxLayers %if at max layers...
    %prob (addLayer) = 0;
    %if numOptions = 4, prob(deleteLayer) = 25%, prob(perturbDepth) = 38%,
    %prob(perturbRho = 38%);
    %if numOptions = 5, prob(deleteLayer) = 20%, prob(perturbDepth) = 27%,
    %prob(perturbRho) = 27%, prob(changeVar) = 27%;
    tmp = rand(); %uniformly distributed
    if tmp < 1/numOptions
        output = 2; %delete layer
    else
        availableOptions = [1 4 5];
        output = availableOptions(randi(length(availableOptions)));
    end
else %If 1 < number of layers < maxLayers...
    output = randi([1,numOptions]); %All options have equal likelihood
end
end