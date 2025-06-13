function [new_depths,new_rhos]=piecewiseLinearSolution(depths,rhos,pbounds)
% interpolate the solution given in depths,rhos onto a uniformly spaced
% vector of depth values.
nlayer = nnz(~isnan(depths));
ndepth = 100;
new_depths = logspace(log10(pbounds.depthMin),log10(pbounds.depthMax),ndepth)';
if nlayer == 1
    new_rhos = ones(ndepth,1)*rhos(1);
else
    new_rhos = rhos(nlayer)*ones(ndepth,1);% assign deep values the last layer value
    log_depth = log(depths(1:nlayer));
    log_depth(1) = log(pbounds.depthMin); % assign top control point location to min depth value

    log_new_depths = log(new_depths);
    mask = log_new_depths < log_depth(nlayer);
    
    lrho = log(rhos(1:nlayer));
    
    new_rhos(mask) = exp(interp1q(log_depth,lrho,log_new_depths(mask)));
    % figure, plot(new_rhos,new_depths)
    % set(gca,'XScale','log','YScale','log');
    % hold on
    % depths(1) = pbounds.depthMin;
    % scatter(rhos,depths)
end