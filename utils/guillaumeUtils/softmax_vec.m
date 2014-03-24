function [res,lss] = softmax_vec(x,inds,lbl)
% softmax_vec - softmax on a vector with non-uniform categories
% 
% [res,lss] = softmax_vec(x,inds)
%
% example:
% softmax_vec([1 2 3 log(1) log(2) 1 1 1 1]',lbl2cell([1 1 1 2 2 3 3 3 3]'))
%
% note: we could use softmax on sparse matrices to do this computation
%   but this vectorial implementation is much facter
% see also lbl2cell_fast, steps


x_cell = lbl2cell_fast(x,inds);
if nargin<3
    lbl = steps(cellfun(@length,inds));
end

maxis = cellfun(@(x) max([max(x),0]),x_cell)';


x_norm_cell = lbl2cell_fast(exp(x-maxis(lbl)),inds);
sums = max(cellfun(@sum,x_norm_cell)',1e-300);

lss = log(sums) + maxis;

res = exp(x-lss(lbl));


