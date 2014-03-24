function [y,grad] = softmax_transfo(x) 
% SOFTMAX_TRANSFO - Softmax bijection with gradient
% 
% [y,grad] = softmax_transfo(x) 
%
% Operates along the columns of x (x has a size (D-1)xN)
% y has a size DxN
% grad is a (D-1)xDxN matrix
% f(x) = exp(x)/(1+sum(exp(x)))
%
% Most useful when maximizing functions of proportions 
% (transform the constraints 0<pi_i<1 and sum(pi_i)==1 into an unconstraint problem.
% 
% softmax_transfo([0 1 -2])
% 
% See also softmax

y = exp(logsoftmax([x;zeros(1,max(1,size(x,2)))]));
[d,n] = size(y);

if nargout>1
    grad = -repmat(reshape(y(1:d-1,:),[d-1,1,n]),[1,d,1]).*...
        repmat(reshape(y(1:d,:)',[1,d,n]),[d-1,1,1]);
        ind = sub2ind([d-1,d,n],(1:d-1)'*ones(1,n),(1:d-1)'*ones(1,n),ones(d-1,1)*(1:n));
    grad(ind) = grad(ind)+y(1:d-1,:);
end

