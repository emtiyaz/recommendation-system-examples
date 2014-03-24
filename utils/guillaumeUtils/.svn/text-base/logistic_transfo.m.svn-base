function [y,grad] = logistic(x) 
% LOGISTIC - Logistic bijection with gradient
% 
% [y,grad] = logistic(x) 
%
% Operates along the columns of x (x has a size 1xN)
% y has a size 1xN
% grad is a 1xDxN matrix
% f(x) = exp(x)/(1+sum(exp(x)))
%
% Most useful when maximizing functions of proportions 
% (transform the constraints 0<pi<1 into an unconstraint problem.
% 
% logistic([0 1 -2])
% 
% See also logistic

y = 1./(1+exp(-x));
[d,n] = size(y);

if nargout>1
    grad = y.*(1-y);
end

