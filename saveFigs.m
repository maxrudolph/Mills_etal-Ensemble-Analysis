function saveFigs(saveFigures,folderName,figNumber)
if saveFigures
    if exist('exportgraphics')
        exportgraphics(gcf,[folderName '/Figure ' figNumber '.png']);
        exportgraphics(gcf,[folderName '/Figure ' figNumber '.eps']);
        exportgraphics(gcf,['Allfigure',figNumber,'/',folderName,'.png']);
        exportgraphics(gcf,['Allfigure',figNumber,'/',folderName,'.eps']);
        close
    else
        export_fig([folderName '/Figure ' figNumber '.png']);
        export_fig([folderName '/Figure ' figNumber '.eps']);
        export_fig(['Allfigure',figNumber,'/',folderName,'.png']);
        export_fig(['Allfigure',figNumber,'/',folderName,'.eps']);
        close
    end
end
end