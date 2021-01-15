function saveFigs(saveFigures,folderName,figNumber)
if saveFigures
    set(gcf,'Position',[100 10 1100 600]);
    %saveas(gcf,[folderName '/Figure ' figNumber],fileType);
    exportgraphics(gcf,[folderName '/Figure ' figNumber '.png']);
    exportgraphics(gcf,[folderName '/Figure ' figNumber '.eps']);
    close
end
end