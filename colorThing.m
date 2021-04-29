exactSlnMisfit = [4.2049, 10.1457,22.5777,37.8140,94.9359;
    4.1909,6.6296,19.8650,32.3442,81.6979;
    3.6012,7.4614,14.5106,20.6446,68.5544];

msMeanMisfit = [3.6867,10.0133,22.5125,37.7310,90.0428;
    19.8255,22.5839,113.2499,91.5695,254.3669;
    18.7288,32.8253,93.9770,100.0531,145.0596];

msMaxLikeMisfit = [3.7591,10.0228,22.4697,37.4412,94.4966;
    19.2358,16.6225,36.9386,40.8314,144.7471;
    9.686,18.3157,29.9407,32.7060,79.0816];

msMedianMisfit = [3.6872,10.0117,22.5136,37.7275,94.2317;
    19.1814,18.1925,34.7947,43.3489,190.9420;
    4.2836,8.8349,45.8229,43.9643,81.9180];

dsMedianMisfit = [3.7264,10.1186,22.7334,37.9893,94.3732;
    3.8030,6.9597,18.0751,32.7155,53.6550;
    3.7377,7.0216,14.6301,21.0194,71.8440];

data1 = msMeanMisfit./exactSlnMisfit;
data1 = (log10(data1)+1)*255/2;
%data1 = data1*40;
data2 = msMaxLikeMisfit./exactSlnMisfit;
%data2 = data2*40;
data2 = (log10(data2)+1)*255/2;

data3 = dsMedianMisfit./exactSlnMisfit;
data3 = (log10(data3)+1)*255/2;
%data3 = data3*40;
tickFontSize = 12
axisFontSize = 14
titleFontSize = 16
cmap = parula

figure
t=tiledlayout(3,1)
nexttile
image(data1)
colormap(cmap);
set(gca,'xticklabels',[])
yticks([1:3])
yticklabels({'1 Layer','3 Layers','4 Layers'});
set(gca,'FontSize',tickFontSize)
title('A: Model Space Mean','FontSize',axisFontSize)

nexttile
image(data2)
colormap(cmap);
set(gca,'xticklabels',[])
yticks([1:3])
yticklabels({'1 Layer','3 Layers','4 Layers'});
set(gca,'FontSize',tickFontSize)
ylabel('Number of layers','FontSize',axisFontSize)
title('B: Maximum Likelihood Model','FontSize',axisFontSize)

nexttile
image(data3)
colormap(cmap);
xticks([1:5])
xticklabels({'0.01','0.02','0.05','0.1','0.2'})
yticks([1:3])
yticklabels({'1 Layer','3 Layers','4 Layers'});
set(gca,'FontSize',tickFontSize)
xlabel('Noise levels','FontSize',axisFontSize)
title('C: Data-Space Median','FontSize',axisFontSize)

cbh = colorbar;
cbh.Layout.Tile = 'east';
cbh.Ticks = linspace(1, 255, 3) ;
cbh.TickLabels = {'0.1','1','10'} ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
 cbh.FontSize = axisFontSize;
 set(gca,'ColorScale','log')
set(get(cbh,'label'),'string','Ratio');
title(t,'Normalized Data Space Misfit Levels','FontSize',titleFontSize)
 
 

%title('Data-Space Misfit Ratio: Maximum Likelihood Model over Exact Solution','FontSize',24)
%xlabel('Noise levels','FontSize',20)
%ylabel('Number of layers','FontSize',20)
%}
%{
data1 = data1*255/max(max(data1));
data2 = data2*255/max(max(data2));
% Display it.
image(data1);
% Initialize a color map array of 256 colors.
colorMap = jet(256);
% Apply the colormap and show the colorbar
colormap(colorMap);
colorbar;

figure
image(data2);
colorMap = jet(256);
colormap(colorMap);
colorbar;
%}