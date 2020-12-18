for i = 1:eva.OptimalK
 centroidYs(:,i) = forwardModel(zVals,10.^GMModel.mu(i,:)',data.lambda);
    % evaluate the forward model for the maximum likelihood.;
    centroidMisfits(i) = norm(centroidYs(:,i)-data.y);
end