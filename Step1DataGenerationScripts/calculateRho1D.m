function apparentResistivity = calculateRho1D(depths,rhos,lambda)
%{
6/16/21 Calculate 1-Dimensional Resisitivity Profile
Simulates Schlumberger measurement of a subsurface structure (assuming
horizontal isotropic homogeneous layers). Takes in subsurface properties
and spacing, outputs the apparent resistivity that would be measured at
those spacings. This fxn uses an 11-point filter from Guptasarma 1982. The
method is from Constable et al 1987 Appendix A.
Inputs:
    depths: column vector of depths to layer interfaces. The first number
        should be 0, representing the surface. Depth is positive in meters.
        Example: [0; 1; 25; NaN; ... NaN;] is a 3-layer model.
    rhos: column vector of resistivity values for each layer, in ohm-meters
    lambda: An 11xm array where m is the number of electrode spacings (so,
        for one measurement, this would be an 11x1 vector). See fxn
        makeLambda. Note the 11-point filter is used in this fxn as well
Output:
    apparentResistivity: a 1xm array of resistivity measurements at each of
    the m electrode spacings.
%}

depths = depths(~isnan(depths));
rhos = rhos(~isnan(rhos));
numLayers = length(depths); %n
h = diff(depths); %layer thicknesses
apparentResistivity = zeros(size(lambda,2),1);

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

%Now calculate T_1(lambda_k) for each k in {1,2,...,11}. Start with
%T_n = rho_n at the top of the terminating half-space (lowest layer
%interface), n being number of layers, then evaluate by Constable et al
%1987 eq A-4 repeatedly to get T_1(lambda_k)

for i = 1:size(lambda,2) %for each electrode spacing...
    T = rhos(numLayers)*ones(size(lambda,1),1); %T_n = rho_n
    for k = (numLayers-1):-1:1 %for each layer above lowest
        T = (T + (rhos(k)*tanh(lambda(:,i)*h(k))))./...
            (1 + (T.*tanh(lambda(:,i)*h(k))/rhos(k))); %Constable A-4
    end
    apparentResistivity(i) = sum(T.*filter); %...the measured resist is:
end
end