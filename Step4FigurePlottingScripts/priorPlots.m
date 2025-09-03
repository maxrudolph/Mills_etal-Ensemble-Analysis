clear;
close all;

filename = 'Ensemble_Constable1984__hierarchical-1_rhoPrior-2_1_18-Dec-2023.mat';
% filename = 'Ensemble_Constable1984_1_PRIOR_04-Dec-2023.mat'
% Ensemble_Constable1984_Wauchope__hierarchical-1_rhoPrior-1_1_PRIOR_30-Nov-2023.mat
load(filename)
%%
nlayer = sum( ~isnan(results.ensembleDepths),1);

figure,
[c,b] = hist(nlayer,1:30);
bar(b,c)
hold on

cc = zeros(size(c));
cc(1) = c(1);
for i=2:30
    cc(i) = cc(1)/(i);
%     cc(i) = cc(i-1)*(i-1)/(i);
end
plot(cc)
%% Histogram of layer depths
mask = ~isnan(results.ensembleDepths);
mask(1,:) = false; % remove the first layer (always at the surface)
logDepths = log10(results.ensembleDepths);

mask = ~isinf(logDepths); % don't include the first layer (depth=0, logDepth = -Inf)
figure, hist(logDepths(mask),100)
xlabel('log10(depth)')


%% look at probability of 2->3 transition and 3->2 transition
mask = results.allChoices==2 & results.allNumLayers == 2;% 2= delete, 3 = add
mask2 = results.allChoices==3 & results.allNumLayers == 1;% 2= delete, 3 = add

figure, 
hist( (results.allProbAccepts(mask)) ,100); title('3->2')
figure, 
hist( (results.allProbAccepts(mask2)),100); title('2->3')


%% compute layer thicknesses
thicknesses = NaN*zeros(size(results.ensembleDepths));
for i = 1:size(results.ensembleDepths,2)
    myDepth = logDepths(:,i);
    myDepth(nlayer(i)+1) = 5;
    myDepth(1) = -1;% surface is effectively at logDepth=-1
    thicknesses(1:nlayer(i),i) = diff(myDepth(1:nlayer(i)+1));%results.ensembleDepths(2:end,i) - results.ensembleDepths(1:end-1,i);
end

% mask = ~isnan(thicknesses) & thicknesses ~= 1e5;
% tmp = thicknesses(mask);
bins = linspace(-1,6.0,100);
figure, hist( thicknesses(:),bins); xlabel('layer thickness')
figure, hist( thicknesses( thicknesses~=6 ),bins); xlabel('layer thickness (no 1-layer models)')

%% histograms of layer thicknesses for each model complexity
figure();
t = tiledlayout(2,5)

for k=2:10
    nexttile();
    m1 = find( nlayer == k );
    % m1 = true(size(nlayer));
    tmp = thicknesses(:,m1);

    hist(tmp(~isnan(tmp)),100); 
    xlabel(['thickness (' num2str(k) '-layer)'])
end

%% histograms of log depths for each model complexity
figure();
t = tiledlayout(2,5)

for k=2:10
    nexttile();
    m1 = find( nlayer == k );
    % m1 = true(size(nlayer));
    tmp = logDepths(2:end,m1);

    hist((tmp(~isnan(tmp))),100); 
    xlabel(['log(depth) (' num2str(k) '-layer)'])
end

%% get the ensemble in long-form
disp('Model space...')
nzplot = 1000; %number of imaginary (depth)layers to divide appraisals into
nRhoBins = 1001; %number of resistivity bins in model space histogram

%Setup logDepthPlot and logRhoPlot
numSavedRuns = size(results.ensembleRhos,2);
minDistL = -1;%log10(min(data.x));
maxDistL = 5;%log10(max(data.x));
zVals = logspace(minDistL,maxDistL,nzplot)'; %depth values for evaluating
logDepthPlot = log10(repmat(zVals,1,numSavedRuns));
logRhoPlot = zeros(length(zVals),numSavedRuns);
for i = 1:numSavedRuns %for each sln...
    logRhoPlot(:,i) = log10(longForm(zVals,results.ensembleDepths(:,i),...
        results.ensembleRhos(:,i)));
end

[numElements,binCenters]=hist3([logRhoPlot(:),logDepthPlot(:)],...
    {linspace(log10(pBounds.rhoMin),log10(pBounds.rhoMax),nRhoBins) ...
    log10(zVals)},'CDataMode','auto');

%%
figure();
figure, surf(binCenters{2},binCenters{1},numElements); shading flat; lighting gouraud;
%%
dx = diff(binCenters{1});
dx = dx(1);
colsum = sum(numElements,2);
figure, plot(binCenters{1},colsum/sum(colsum)/dx);
hold on
plot(binCenters{1},normpdf(binCenters{1},3,1),'r--');
