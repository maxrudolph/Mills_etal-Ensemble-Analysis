function output = chooseOption(currentLayers,maxLayers,numOptions)
if currentLayers == 1 %If one layer model...
    % prob(remove layer) = prob(perturbDepth) = 0;
    %numOptions = [4,5] => prob(addLayer) = [25%,20%];
    %prob(perturbRho) = [75%,48%],prob(changeVar) = [0%,32%]
    tmp = rand(); %uniformly distributed
    if tmp<(1/numOptions)
        output = 3; %Add new layer.
    elseif (tmp < 1/numOptions + (1-1/numOptions)/(numOptions-3))
        output = 4; %perturb resistivity
    else
        output = 5; %var
    end
elseif currentLayers >= maxLayers %If at max layers...
    %prob (addLayer) = 0; numOptions = [4,5] =>
    %prob(deleteLayer) = [25%,20%]; prob(perturbDepth) = [46.9%,37.3%];
    %prob(perturbRho) = [28.1%,31.3%], prob(changeVar) = [0%,11.4%]
    tmp = rand(); %uniformly distributed
    if tmp < 1/numOptions
        output = 2; %delete layer
    elseif (tmp < 1/numOptions + (1-1/numOptions)/(numOptions-2))
        output = 1; %perturb depth
    elseif (tmp < 1/numOptions + (2-2/numOptions)/(numOptions-2))
        output = 4; %perturb resistivity
    else
        output = 5;
    end
else %normal circumstances
    output = randi([1,numOptions]);
end
end