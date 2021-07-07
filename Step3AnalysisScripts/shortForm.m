function [outDepths,outRhos] = shortForm(inDepths,inRhos)
%removes duplicate values
h = diff(inRhos);
ind = find(h==0);
outRhos = inRhos;
outDepths = inDepths;
outRhos(ind+1) = [];
outDepths(ind+1) = [];
end