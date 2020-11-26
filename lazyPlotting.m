function lazyPlotting(subplotSize1,subplotSize2,histTF,things,titles,...
    xlabels,ylabels)
    figure();
    for i = 1:(subplotSize1*subplotSize2)
        subplot(subplotSize1,subplotSize2,i);
        if histTF(i)
            histogram(things{i})
        else
            plot(things{i})
        end
        title(titles{i});
        xlabel(xlabels{i});
        ylabel(ylabels{i});
    end
end
        