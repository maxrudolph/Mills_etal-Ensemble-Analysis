function filename = loadFieldDataRenner(varargin)
%{
Saves in a file these things:
    forwardModel: Function handle of the forward model
    data: a structure containing data and related information
        x: values of electrode spacings, size numMeasurements, log-spaced
            from minDist to maxDist, in meters.
        lambda: a 2D array based on x, see makeLambda and calculateRho1D
        fx: the measurements with NO noise
        y: the measurements WITH noise
        Cd: covariance matrix, which in our case is
Then outputs a filename which can easily be fed into the inversion script,
which is the next step of the process.
%}

%% 0 Preliminary stuff ignore this section
addpath(genpath(fileparts(mfilename('fullpath'))))
%adds subfolders so you can use the scripts in them
defaultNoise = 1.0;
defaultSubStruct = 'Constable1984_Renner';

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x>=0);
addOptional(p,'noise',defaultNoise,validScalarPosNum);
addParameter(p,'subStructChoice',defaultSubStruct,@ischar);
parse(p,varargin{:});

%% 1 User Set Options Here:

filterSize = 11;
%Choice of filter size for the forward model. Choices are 7,11,19. Based on
%Guptasarma 1982. This effects both which script is used for the forward
%model as well as the size of the lambda matrices, so it is an accuracy vs.
%computational expense tradeoff.

data.subStructChoice = 'Constable1984_Renner'; 
%see subStructGen for choices

% minDist = 0.1; maxDist = 1000; %in meters
%numMeasurements = 21; %how many measurements
% data.x = logspace(log10(minDist),log10(maxDist),numMeasurements);
% array of electrode spacings

data.noiseCoef = 1.0; %How "noisy" are the measurements.


%% Digitized data from Constable paper
% digitization was done with web plot digitizer.

% This is a 1-sigma standard deviation.
error_bar = [10.131384002971298, 8.23443115885666
14.144866165713562, 7.857929039012548
20.24396481119087, 7.777796778179991
28.272693327019812, 8.904672537539831
40.473215678047694, 10.073107632971464
55.846095158455526, 13.503255184158245
80.98209114229194, 20.443313299857543
101.22527627690566, 28.4064648320283
150.4575778484298, 42.49268347230918
202.56614513745888, 62.00686067997087
6093.604993572734, 626.2086253319778
7615.175569638835, 770.6556653769237
10246.556557335902, 810.2697195608571
14314.369395778906, 1086.257825836283
20241.875620838753, 1305.6068454286926
27583.543543943226, 1627.1379173244823
40000.470973162955, 2523.9578034523297
60925.12797621269, 3185.820268804728
92829.97002504185, 4943.018425342695
];

apparent_resistivity_data = [10.12742243056727, 6.618039705423816
13.966213791244554, 6.470241253607495
20.237368102557017, 6.482906070623123
28.263480372842473, 7.422173309089261
39.96203119609014, 8.294225645879676
55.825471164753864, 10.985167865085513
80.95570227702724, 17.039797217009504
99.94678160956249, 23.389964418466732
150.40528168845438, 34.99087787630361
202.49573704181398, 51.05995461453317
6091.354621892286, 509.4332273233291
7612.197884716454, 619.3785210273652
10242.772498355554, 659.1706047156835
14309.393996316985, 894.4861049028713
20235.27959298953, 1088.2421829048337
27916.96852862111, 1339.9646903523553
39986.56757887064, 2078.369546208952
60886.75116641176, 2240.371317573986
91597.78429448836, 2827.6840136107853
    ];

figure()
plot(apparent_resistivity_data(:,1),apparent_resistivity_data(:,2),'o')
hold on
plot(apparent_resistivity_data(:,1),error_bar(:,2),'r.')
set(gca,'XScale','log','YScale','log')

relative_error = (error_bar(:,2)-apparent_resistivity_data(:,2))./apparent_resistivity_data(:,2);
% note - all but the last TWO data points are close to 20% relative error.
% The paper states that the relative error is at least 20% and this appears
% to be the value adopted for most of the data.
relative_error(1:end-2) = 0.2;
data.relative_error = relative_error;
%% 2 Calculated Stuff:

switch filterSize
    case 7
        forwardModel = @(a,b,c) calculateRho1D07(a,b,c);
    case 11
        forwardModel = @(a,b,c) calculateRho1D11(a,b,c);
    case 19
        forwardModel = @(a,b,c) calculateRho1D19(a,b,c);
end %anything else is an invalid choice

% [trueDepths,trueRhos] = subStructGen(data.subStructChoice);
data.x = apparent_resistivity_data(:,1)';
data.lambda = makeLambda(data.x,filterSize); %lambda matrix for calculateRho1D
% data.fx = forwardModel(trueDepths,trueRhos,data.lambda); %'true' output
data.fx = apparent_resistivity_data(:,2);
data.y = apparent_resistivity_data(:,2);

rng(1); %re-seed random number generator for consistent noise pattern.
% noiseVector = data.noiseCoef.*data.fx.*randn(length(data.fx),1);
%noise is added, Gaussian with mean 0 and std dev = noiseCoef*f(x)

data.Cd = diag( (relative_error .* apparent_resistivity_data(:,2)).^2 );

%% 3 Save
filename = ['Data_', data.subStructChoice, '_',...
    num2str(data.noiseCoef), '_', date, '.mat'];
save(filename,'data','forwardModel');
end