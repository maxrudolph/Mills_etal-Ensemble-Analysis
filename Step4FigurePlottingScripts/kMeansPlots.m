function kMeansPlots(intitle,idx,sumd,inModels,visibility)
numClusters = size(inModels,2)-1;
figure('visible',visibility,'units','normalized','outerposition',[0 0 1 1]);
title(intitle)
subplot(3,1,1);
histogram(idx,'EdgeAlpha',0);
%set(gca,'YScale','log');
title('Number of ensemble slns');
xlabel('Cluster #');
ylabel('Frequency');

averageDist = zeros(1,numClusters);
centroidMisfits = zeros(1,numClusters);
for i = 1:numClusters
    averageDist(i) = sumd(i)/nnz(idx==i);
    centroidMisfits(i) = inModels{i+1}.misfit;
end

subplot(3,1,2);
bar(averageDist);
title('Distance from Centroid ');
ylabel('Log-distance');
xlabel('Cluster');

subplot(3,1,3);
bar(centroidMisfits)
title('Misfit')
ylabel('Data-space misfit (m)')
xlabel('Cluster');
end