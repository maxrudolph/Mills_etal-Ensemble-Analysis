function saveFigs(saveFigures,folderName,figNumber)
if saveFigures
    disp(['Saving figure ' num2str(figNumber)])
    if exist('exportgraphics')
        dpi = get(groot,'ScreenPixelsPerInch');  
        set(gcf,'Units','pixels','Position',[0 0 7*dpi 4.5*dpi]);
        exportgraphics(gcf,[folderName '/Figure ' sprintf(figNumber) '.png']);
        %exportgraphics(gcf,[folderName '/Figure ' sprintf(figNumber) '.eps'],'ContentType','vector');
        exportgraphics(gcf,['Allfigure',figNumber,'/',folderName,'.png']);
        %exportgraphics(gcf,['Allfigure',figNumber,'/',folderName,'.eps'],'ContentType','vector');
        close
    else
        export_fig([folderName '/Figure ' figNumber '.png']);
        %export_fig([folderName '/Figure ' figNumber '.eps']);
        export_fig(['Allfigure',figNumber,'/',folderName,'.png']);
        %export_fig(['Allfigure',figNumber,'/',folderName,'.eps']);
        close
    end
end
end