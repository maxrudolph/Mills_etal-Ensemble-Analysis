function lambda = makeLambda(x)
%{
6/16/21 Make Lambda Matrix
Input: x, an array of n electrode spacings as in a Schlumberger array
Output: an 11xn vector of 'lambda' values, based on an 11-point filter from
Guptasarma 1982. These lambda values are then used in calculateRho1D. This
fxn is defined separately because it only needs to run once for a given
electrode spacing, so this cuts down on overall runtime
%}

%See Guptasarma 1982 pg 512 for abciss values
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

[xGrid,abcissGrid] = meshgrid(x,abciss);
lambda = 10.^(abcissGrid-log10(xGrid)); %Guptasarma 1982 pg 512
end