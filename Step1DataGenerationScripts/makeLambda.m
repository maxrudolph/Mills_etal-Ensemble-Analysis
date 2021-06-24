function lambda = makeLambda(x,filterSize)
%{
6/16/21 Make Lambda Matrix
Input: x, an array of n electrode spacings as in a Schlumberger array
filterSize, choice of filter based on Guptasarma 1982, should match the
forward model script being used as well.
Output: a 7,11,or 19xn vector of 'lambda' values, based on the filters from
Guptasarma 1982. These lambda values are then used in calculateRho1D. This
fxn is defined separately because it only needs to run once for a given
electrode spacing, so this cuts down on overall runtime
%}

%See Guptasarma 1982 pgs 512,513 for these values
switch filterSize
    case 7
        abciss = [
            -0.17445;
            0.09672;
            0.36789;
            0.63906;
            0.91023;
            1.1814;
            1.45257];
    case 11
        abciss = [
            -0.420625;
            -0.20265625;
            0.0153125;
            0.23328125;
            0.45125;
            0.66921875;
            0.8871875;
            1.10515625;
            1.323125;
            1.54109375;
            1.7590625];
    case 19
        abciss = [
            -0.980685;
            -0.771995;
            -0.563305;
            -0.354615;
            -0.145925;
            0.062765;
            0.271455;
            0.480145;
            0.688835;
            0.897525;
            1.106215;
            1.314905;
            1.523595;
            1.732285;
            1.940975;
            2.149665;
            2.358355;
            2.567045;
            2.775735];
end

[xGrid,abcissGrid] = meshgrid(x,abciss);
lambda = 10.^(abcissGrid-log10(xGrid));
end