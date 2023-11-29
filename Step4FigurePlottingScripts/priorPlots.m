clear;
% close all;

filename = 'Ensemble_Constable1984_Wauchope__hierarchical-1_rhoPrior-2_1_PRIOR_27-Nov-2023.mat';
load(filename)

nlayer = sum( ~isnan(results.ensembleDepths),1);
figure, hist(nlayer,1:30);

for i=1:length(results.ensembleDepths)
   results.ensembleDepths(nlayer(i)+1,i) = 1e5; 
end
mask = ~isnan(results.ensembleDepths);
mask(1,:) = false;

%% hist on layer depths
allDepths = results.ensembleDepths;
mydata = log10(allDepths);
mask1 = ~isinf(mydata) & mydata ~= 5; % don't include bottom layer position.
mydata = mydata(mask1);
figure, hist(mydata(:),100);

%% compute layer thicknesses
thicknesses = NaN*zeros(size(results.ensembleDepths)-[1 0]);
for i = 1:length(results.ensembleDepths)
    thicknesses(:,i) = results.ensembleDepths(2:end,i) - results.ensembleDepths(1:end-1,i);
end

mask = ~isnan(thicknesses) & thicknesses ~= 1e5;
tmp = thicknesses(mask);
figure, hist(log10(tmp(:)),100);

%% 1-layer models
m1 = nlayer == 1;
tmp = thicknesses(:,m1);

figure, hist(log10(tmp(~isnan(tmp))),100)