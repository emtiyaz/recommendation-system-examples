function x = logistic_transfo_inv(y) 
% LOGISTIC_TRANSFO_INV - Inverse logistic (log-odds)
% 
% [y,grad] = logistic(x) 
%
% Operates along the columns of x (x has a size 1xN)
% y has a size 1xN
% grad is a 1xDxN matrix
% f(x) = log(x/(1-x));
%
% Most useful when maximizing functions of proportions 
% (transform the constraints 0<pi<1 into an unconstraint problem.
% 
% logistic_inv([0 .1 .5 1])
% 
% See also logistic_transfo

x = log(max(1e-320,y./max(1e-320,(1-y))));

