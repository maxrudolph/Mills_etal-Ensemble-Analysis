function [depths,rhos] = modelGen(kMax,choice)
    if choice == '3LayerA'
        inDepths = [0,1,25]';
        inRhos = [10,390,10]';
    elseif choice == '4LayerA'
        inDepths = [0,5,30,200]';
        inRhos = [100,10,250,1]';
    elseif choice == '4LayerB'
        inDepths = [0,5,30,200]';
        inRhos = [100,10,250,1e8]';
    elseif choice == '4LayerC'
        inDepths = [0,5,30,200]';
        inRhos = [100,10,250,1e5]';
    elseif choice == '4LayerD'
        inDepths = [0,5,30,200]';
        inRhos = [100,10,250,1e2]';
    elseif choice == '1LayerA'
        inDepths = [0]';
        inRhos = [100]';
    end
    
    depths = nan*zeros(kMax,1);
    rhos = nan*zeros(kMax,1);
    depths(1:size(inDepths)) = inDepths;
    rhos(1:size(inDepths)) = inRhos;
end
