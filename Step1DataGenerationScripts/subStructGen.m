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

    if choice == '3LayerA'
        depths = [0,1,25]';
        rhos = [10,390,10]';
    elseif choice == '4LayerA'
        depths = [0,5,30,200]';
        rhos = [100,10,250,1]';
    elseif choice == '4LayerB'
        depths = [0,5,30,200]';
        rhos = [100,10,250,1e-2]';
    elseif choice == '4LayerC'
        depths = [0,5,30,200]';
        rhos = [100,10,250,1e2]';
    elseif choice == '4LayerD'
        depths = [0,5,30,200]';
        rhos = [100,10,250,1e4]';
    elseif choice == '1LayerA'
        depths = [0]';
        rhos = [100]';
    end
end
