function [res,lss] = logsoftmax(x,dim) 
% LOGSOFTMAX - Logarithm of the Softmax activation function
% 
% s = logsoftmax(x) 
% where the rows of x are normalized
%
% softmax = exp(x)/(sum(exp(x),2))
% 
% [res,ls] = logsoftmax(x) also returns the log sum of the exponential of each row
%
% see also normalizedim, softmax, softmax_vec

if nargin<2
    dim=1;
end

sz=size(x);
rep = ones(1,length(sz));
rep(dim) = sz(dim);
maxi = max(x,[],dim);
lss = log(sum(exp(x-repmat(maxi,rep)),dim)) + maxi;
res = x-repmat(lss,rep);


