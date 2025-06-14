function apparentResistivity = piecewiseLinearWrapper(depths,rhos,lambda,forwardModel,pBounds)

[new_depths,new_rhos]=piecewiseLinearSolution(depths,rhos,pBounds);
apparentResistivity = forwardModel(new_depths,new_rhos,lambda);

