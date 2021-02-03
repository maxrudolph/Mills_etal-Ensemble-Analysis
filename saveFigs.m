function saveFigs(saveFigures,folderName,figNumber)
if saveFigures
    exportgraphics(gcf,[folderName '/Figure ' figNumber '.png']);
    exportgraphics(gcf,[folderName '/Figure ' figNumber '.eps']);
    exportgraphics(gcf,['Allfigure',figNumber,'/',folderName,'.png']);
    exportgraphics(gcf,['Allfigure',figNumber,'/',folderName,'.eps']);
    close
end
end