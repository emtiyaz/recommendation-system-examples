function [nlZ] = computeMargLikLowerBoundKL(post, hyp, mean, cov, lik, x, y)
% compute KL lower bound given the posterior distribution
% post is posterior distribution, rest is same input as other gpml functions
% Also, y is -1 or 1, 
% For usage, see testInferBernLogitGP_pw.m
% Written by Emtiyaz, CS, UBC, 
% Modified on Sep 27, 2011

  alpha = post.alpha;
  lambda = (post.sW).^2;

  y = (y+1)/2; % convert to 0,1 from -1 +1
  covMat = feval(cov{:}, hyp.cov, x); 
  mean_ = feval(mean{:}, hyp.mean, x);
  bound = getbound(20);

  nlZ=funObj_BernLogitGP_pw([alpha(:); log(lambda(:))],y,mean_,covMat,bound);

