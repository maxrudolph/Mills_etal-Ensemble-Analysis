function apparentResistivity = calculateRho1D19(depths,rhos,lambda)
%{
6/16/21 Calculate 1-Dimensional Resisitivity Profile w/19 point filter
Simulates Schlumberger measurement of a subsurface structure (assuming
horizontal isotropic homogeneous layers). Takes in subsurface properties
and spacing, outputs the apparent resistivity that would be measured at
those spacings. This fxn uses a 19-point filter from Guptasarma 1982. The
method is from Constable et al 1987 Appendix A.
Inputs:
    depths: column vector of depths to layer interfaces. The first number
        should be 0, representing the surface. Depth is positive in meters.
        Example: [0; 1; 25; NaN; ... NaN;] is a 3-layer model.
    rhos: column vector of resistivity values for each layer, in ohm-meters
    lambda: A 19xm array where m is the number of electrode spacings (so,
        for one measurement, this would be an 11x1 vector). See fxn
        makeLambda. Note the 19-point filter is used in this fxn as well
Output:
    apparentResistivity: a 1xm array of resistivity measurements at each of
    the m electrode spacings.
%}

depths = depths(~isnan(depths));
rhos = rhos(~isnan(rhos));
numLayers = length(depths); %n
h = diff(depths); %layer thicknesses
apparentResistivity = zeros(size(lambda,2),1);

%filter coefficients (CONSTANTS)19-point filter, Guptasarma 1982 table A3
filter = [
    0.00097112;
   -0.00102152;
    0.00906965;
    0.01404316;
    0.09012;
    0.30171582;
    0.99627084;
    1.3690832;
   -2.99681171;
    1.65463068;
   -0.59399277;
    0.22329813;
   -0.10119309;
    0.05186135;
   -0.02748647;
    0.01384932;
   -0.00599074;
    0.00190463;
   -0.0003216]; 
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