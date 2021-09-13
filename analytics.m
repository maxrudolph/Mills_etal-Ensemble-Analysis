filenames = {
    '3LayerA_0_02-Jul-2021.mat';
    '3LayerA_0.01_02-Jul-2021.mat';
    '3LayerA_0.02_02-Jul-2021.mat';
    '3LayerA_0.05_02-Jul-2021.mat';
    '3LayerA_0.1_02-Jul-2021.mat';
    '3LayerA_0.2_02-Jul-2021.mat'};
titles = {'0','0.01','0.02','0.05','0.1','0.2'};
numEnsembles = length(filenames);

for i = 1:numEnsembles    
    load(['Ensemble_' filenames{i}],...
        'results','data','forwardModel');
    titles{i}
    y1(i) = sqrt(10.^mean(log10(results.ensembleVars)));
    y2(i) = sqrt(std(log10(results.ensembleVars)));
end
x = [0,0.01,0.02,0.05,0.1,0.2];
figure, plot(x,y2)
figure, plot(x, (abs(y1-x)./abs(x)))



%{

1.0695e-4+1.6614e-9
.0111+4.9812e-5
.0195e-4+1.5828e-4
.0415+7.063e-4
.1294+.0067
.2457+.0256

%}