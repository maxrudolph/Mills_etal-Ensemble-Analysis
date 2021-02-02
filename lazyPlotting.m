function lazyPlotting(subplotSize1,subplotSize2,histTF,things,titles,...
    xlabels,ylabels,visibility)
    figure('visible',visibility,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:(subplotSize1*subplotSize2)
        subplot(subplotSize1,subplotSize2,i);
        if histTF(i)
            histogram(things{i},'EdgeAlpha',0)
        else
            plot(things{i})
        end
        title(titles{i});
        xlabel(xlabels{i});
        ylabel(ylabels{i});
    end
end
        