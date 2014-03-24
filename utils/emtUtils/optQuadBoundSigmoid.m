function [xi, A, b, c] = optQuadBoundSigmoid(mean_, covMat)
% [xi, A, b, c] = optQuadBoundSigmoid(mean_, covMat)
% Optimize the upper bound to E(log(1+exp(z))). Note that the upper bound is on
% log(1+exp(z)), not the usual definition of sigmoid. 
% (mean_ , covMat) define the distribution wrt which expectation is taken. If
% mean_ is vector (and covMat is matrix of corresponding dimension) then the
% result is the optimization over each dimension separately.
% 
% Written by Emtiyaz, CS, UBC.
% Modified on Feb. 1, 2010

  mean_ = mean_(:);
  % optimize variational parameter
  xi = sqrt((mean_).^2 + diag(covMat));
  [A, b, c] = computeQuadBoundSigmoid(xi);
 
