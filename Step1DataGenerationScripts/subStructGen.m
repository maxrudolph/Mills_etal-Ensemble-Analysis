function [depths,rhos] = subStructGen(choice)
%{
6/16/21 Subsurface Structure Generator
Basically just a database of subsurface structures, for ease of synthetic
data generation pre-inversion. To add a new one, give it a column vector of
positive depths (in meters) to layer interfaces (starting with 0 for top 
layer) and a column vector of resistivities in ohm-meters for each layer.
Title it by its number of layers and A,B,C...

Inputs:
    choice: the string representing which model you're using. 

Outputs:
    depths and rhos are the depth and rho vectors associated with choice,
    formatted appropriately for the other functions that will use them.
%}

    if strcmp(choice,'3LayerA')
        depths = [0,1,25]';
        rhos = [10,390,10]';
    elseif strcmp(choice,'4LayerA')
        depths = [0,5,30,200]';
        rhos = [100,10,250,1]';
    elseif strcmp(choice,'4LayerB')
        depths = [0,5,30,200]';
        rhos = [100,10,250,1e-2]';
    elseif strcmp(choice,'4LayerC')
        depths = [0,5,30,200]';
        rhos = [100,10,250,1e2]';
    elseif strcmp(choice,'4LayerD')
        depths = [0,5,30,200]';
        rhos = [100,10,250,1e4]';
    elseif strcmp(choice,'1LayerA')
        depths = [0]';
        rhos = [100]';
    elseif strcmp(choice,'8LayerA')
        depths = [0,2,7,8,18,35,39,109]';
        rhos = [48,51,923,6220,2,5,19,10]';
    elseif strcmp(choice,'8LayerB')
        depths = [0,1.6,2.4,3.9,4.1,6.6,18.8,62.6]';
        rhos = [16,3,14,1190,366,1604,117,9]';
    elseif strcmp(choice,'6LayerA')
        depths = [0,2,8,35,83,202]';
        rhos = [112,168,888,17,39,109]';
    else
        depths = [];
        rhos = [];
    end
end
