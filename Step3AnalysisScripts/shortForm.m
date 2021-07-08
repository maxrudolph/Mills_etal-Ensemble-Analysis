function [outDepths,outRhos] = shortForm(inDepths,inRhos)
%{
Short Form 7/7/21
Takes a solution of 'long form' (ie see longForm script) and deletes the
depths (and corresponding resistivities) which have the same resistivity as
the depth value just above it. Thus it only keeps the depths which indicate
layer interfaces. This makes it the appropriate format to be put through
the forward model (see any version of calculateRho1D script to see why)
%}
h = diff(inRhos);
ind = find(h==0)+1; 
%if element h_i = 0 (ith element of h), then there is no resistivity
%difference between inRhos_i and inRhos_{i+1}. Similarly, h_i != 0
%indicates inRhos_{i+1} is the start of a new resistivity value and thus
%inDepths_{i+1} is where a new layer starts. (technically the interface 
%could exist anywhere between inDepths_i and inDepths_{i+1}). So if we find
%all the places h=0, then those are the redundant values that can be
%deleted.
outRhos = inRhos;
outDepths = inDepths;
outRhos(ind) = [];
outDepths(ind) = [];
end