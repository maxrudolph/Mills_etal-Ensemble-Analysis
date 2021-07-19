function figurePlotting(filename, saveFigures)
addpath(genpath(fileparts(mfilename('fullpath'))))
filename1 = ['Ensemble_' filename];
filename2 = ['Analysis_' filename];
load(filename1,'data','results')
load(filename2,'allModels','binCenters','KMModelsEuclid','KMModelsMan',...
    'numElements','xVals','yVals');

if saveFigures
    slashpos = find(filename == '/',1,'last');
    ensembleName = filename(slashpos+10:end-9);
    folderName = ['figures_' ensembleName];
    mkdir(folderName);
else
    folderName = ' ';
end

%% Figures 1 through 3
disp('Producing first three figures...')
firstThreeFigures(results,saveFigures,folderName)

%% Figure 4
disp('Figure 4...')
bigPlot(binCenters,numElements,allModels,xVals,yVals,data,results,' ',...
    saveFigures,folderName,'4');

%% Figures 5 and 6
disp('Cluster plotting...')

bigPlot(binCenters,numElements,KMModelsEuclid,xVals,yVals,data,results,...
    'K-means: Euclidean',saveFigures,folderName,'5');
bigPlot(binCenters,numElements,KMModelsMan,xVals,yVals,data,results,...
    'K-means: Manhattan',saveFigures,folderName,'6');
end