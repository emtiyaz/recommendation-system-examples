function [A, b, c] = computeQuadBoundSigmoid(xi)

  lambda = (1./(1+exp(-xi))-0.5)./(2*max(xi,1e-10));
  % compute quadBound
  A = 2*diag(lambda);
  b = -0.5*ones(length(lambda),1);
  c = sum(-lambda.*(xi.^2) - 0.5*xi + log(1+exp(xi)));
 
