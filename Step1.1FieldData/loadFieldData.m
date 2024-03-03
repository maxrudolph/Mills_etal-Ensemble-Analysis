function filename = loadFieldData(varargin)
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
defaultSubStruct = 'Constable1984_Wauchope';

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

data.subStructChoice = p.Results.subStructChoice;%'';%'Constable1984_Wauchope';

% minDist = 0.1; maxDist = 1000; %in meters
%numMeasurements = 21; %how many measurements
% data.x = logspace(log10(minDist),log10(maxDist),numMeasurements);
% array of electrode spacings

data.noiseCoef = 1.0; %How "noisy" are the measurements. (1.0 for field - actual uncertainties are given).


%% Digitized data from Constable paper
% digitization was done with web plot digitizer.

% Wauchope:
if strcmp(data.subStructChoice,'Constable1984_Wauchope')
    % digitization was done with web plot digitizer.

    % This is a 1-sigma standard deviation.
    error_bar = [4.896732901846811, 934.0505282049073
        6.801471165208465, 384.755808736745
        9.945979346038573, 116.59144011798311
        13.807062345759334, 38.47558087367446
        19.879526249467705, 24.694065319965773
        27.65624883320632, 19.119717955005
        39.848815103186915, 16.39890367222246
        54.52034719114441, 16.967959918688983
        80.00815620613719, 20.821811885006582
        100.36068813004117, 23.46228848142265
        152.49729953505766, 28.791166380223558
        201.5900917326827, 35.33036694992737
        295.92082109390674, 48.852735715193866
        499.59268447829635, 84.31909292866251
        986.8252174635064, 163.98903672222477
        1474.4160092126976, 251.18864315095797
        1949.5711452737073, 341.4548873833601
        2912.8551407301143, 523.0202675339294
        3919.72307148123, 748.2971205750629
        5656.767398881645, 1206.3726468100235
        7610.142656820342, 1558.0901878706047
        9887.254908078085, 1978.3188827841643
        33498.83900757835, 3593.813663804626
        50031.26946400382, 4721.435640802164
        69621.08503431256, 4049.5559845768885
        98642.22776801953, 4410.059454176737];

    apparent_resistivity_data = [4.895470008598941, 843.1909292866251
        6.799717029914117, 347.3287559085802
        9.77136924386588, 103.4700871341198
        13.803501431541855, 34.73287559085806
        19.874399214774463, 22.291954510236945
        27.649116140733, 17.259850793256216
        39.83853788355126, 14.803703235666639
        54.506286105984294, 15.317404637020783
        79.99096039856323, 19.119717955005
        100.33049130302543, 20.821811885006582
        152.4579696770154, 25.99051079393097
        201.53810058454602, 31.89361179185746
        295.84450155303, 44.100594541767364
        499.442365522213, 74.82971205750628
        986.5707100339717, 148.0370323566664
        1499.9246871199446, 226.7543125870802
        1949.0683405904724, 308.23992397451434
        2911.978712076172, 464.15888336127773
        3918.5436933627207, 664.0827850634838
        5655.065373618616, 1070.6058931542352
        7608.1799608822075, 1406.5272421052364
        9714.093360746945, 1785.878265500129
        33488.759775703926, 3189.361179185746
        49152.9280992172, 4190.079105786669
        69603.1293891994, 3655.6361467894103
        100314.31816901018, 3473.2875590858057];

    figure()
    plot(apparent_resistivity_data(:,1),apparent_resistivity_data(:,2),'o')
    hold on
    plot(apparent_resistivity_data(:,1),error_bar(:,2),'r.')
    set(gca,'XScale','log','YScale','log')

    relative_error = (error_bar(:,2)-apparent_resistivity_data(:,2))./apparent_resistivity_data(:,2);
    % note - all but the last data point are close to 10% relative error.
    % The paper states that the relative error is at least 10% and this appears
    % to be the value adopted for most of the data.
    relative_error(1:end-1) = 0.1;
    data.relative_error = relative_error;
elseif strcmp(data.subStructChoice,'Constable1984_Renner')
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
else
    error('undefined dataset')
end



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
