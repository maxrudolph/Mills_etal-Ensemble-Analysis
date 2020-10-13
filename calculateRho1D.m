function apparentResistivity = calculateRho1D(depths,rhos,lambda)
    %{
    depths (column vector; depths to layer interfaces)
    rhos (column vector; resistivities of each layer), 
    Each column vector has length = max #oflayers; if numLayers < maxLayers
    , NaN values for missing layers;
    lambda is a meshgrid based on the x values and abcissae. It is defined
    outside this function to cut down on runtime. The lambda values are
    from the first equation, page 512, Guptasarma 1982. Note this fxn
    relies on using an 11-point filter
    %}

    depths = depths(~isnan(depths));
    rhos = rhos(~isnan(rhos));
    numLayers = length(depths);
    h = diff(depths); %layer thicknesses

    apparentResistivity = zeros(size(lambda,2),1); 
    %pre-allocating 'measured' apparent resistivities

%filter coefficients (CONSTANTS)11-point filter, Guptasarma 1982 table a-2
   filter = [
        0.041873;
        -0.022258;
        0.387660;
        0.647103;
        1.84873;
        -2.96084;
        1.358412;
        -0.377590;
        0.097107;
        -0.024243;
        0.004046];

    %Now figure out T_1(lambda_k) for each k. Start with T_n = rho_n at the
    %top of the terminating half-space, then evaluate Constable et al 1987
    %eq A-4 repeatedly to get T_1(lambda_k)

    for i = 1:size(lambda,2)

        T = rhos(numLayers)*ones(size(lambda,1),1);
        
        for k = (numLayers-1):-1:1
            T = (T + (rhos(k)*tanh(lambda(:,i)*h(k))))./...
                (1 + (T.*tanh(lambda(:,i)*h(k))/rhos(k)));
        end
        apparentResistivity(i) = sum(T.*filter);
    end
end