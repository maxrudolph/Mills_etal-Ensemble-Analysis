function lambda = makeLambda(x)
    %See Guptasarma 1982 for abciss values
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
    lambda = 10.^(abcissGrid-log10(xGrid));
 end