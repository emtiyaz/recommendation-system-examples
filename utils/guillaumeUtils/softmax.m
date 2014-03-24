function [sm,s] = softmax(x,dim) 
% SOFTMAX - Softmax activation function
% 
% y=softmax(x) 
%
% Return the exp(x)/sum(exp(x))
%
% [y,s] = softmax(x)
% also return the log-sum in the denominator (log of the normalizing constant)
% 
% See also logodds, logsoftmax
if nargin<2
    dim=1;
end

[sm,s] = logsoftmax(x,dim);
sm = exp(sm);
