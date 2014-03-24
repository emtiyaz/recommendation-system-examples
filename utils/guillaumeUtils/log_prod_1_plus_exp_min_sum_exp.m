function [f,grad] = log_prod_1_plus_exp_min_sum_exp(x,withone)
% log_prod_1_plus_exp_min_sum_exp - log(prod(1+x_k)-sum(exp(x_k))-1)
%  
% usage:
% [f,grad] = log_prod_1_plus_exp_min_sum_exp(x)
%
% example 1:
% checkderivatives(@(x) log_prod_1_plus_exp_min_sum_exp(x),randn(6,1))
%
% example 2:
% pt=-10:.5:10;[x,y]=meshgrid(pt,pt);z=log_prod_1_plus_exp_min_sum_exp([x(:),y(:)]',0);surf(x,y,reshape(z,size(x)),'EdgeColor','none')
% 
if nargin<2
    withone = 1;
end
if withone
    x = [0;x];
end
    
a = sum(logoneplusexp(x));
b = logsumexp(x,1);


v = log_exp_min_one(a-b);

f = b + v;

if nargout>1
    l1pe = logoneplusexp(-x);%log(da/dx)
    lsm = logsoftmax(x,1);%
    loggrad = lsm-v + log_exp_min_one(a-b-l1pe-lsm);
    grad = exp(loggrad);
    if withone
        grad = grad(2:end);
    end
end

