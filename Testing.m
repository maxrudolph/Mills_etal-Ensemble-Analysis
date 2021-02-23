measure.modelChoice = '4LayerA'; %currently setup: 3LayerA, 4LayerA
    %Measurement options
    measure.minDist = 0.1; % Smallest electrode distance, meters
    measure.maxDist = 1000; %  Largest electrode distance, meters
    measure.numMeasurements = 21; %total # of measurements
    measure.noiseCoef = 0.01; %How "noisy" are the measurements
    measure.kMax = 10;
    forwardModel = @(a,b,c) calculateRho1D(a,b,c);
    
    data1 = createSyntheticData(measure,forwardModel);
    measure.modelChoice = '4LayerD';
   data2 = createSyntheticData(measure,forwardModel);
   measure.modelChoice = '4LayerC';
   data3 = createSyntheticData(measure,forwardModel);
   measure.modelChoice = '4LayerB';
   data4 = createSyntheticData(measure,forwardModel);
   
   figure,plot(data1.x,data1.y);
   set(gca,'XScale','log','YScale','log');
   hold on
   plot(data2.x,data2.y);
   set(gca,'XScale','log','YScale','log');
   plot(data3.x,data3.y);
   set(gca,'XScale','log','YScale','log');
   plot(data4.x,data4.y);
   set(gca,'XScale','log','YScale','log');
  
   %....