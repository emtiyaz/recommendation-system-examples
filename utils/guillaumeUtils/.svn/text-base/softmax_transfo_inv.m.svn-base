function x = softmax_transfo_inv(p) 
% SOFTMAX_TRANSFO_INV - inverse softmax transformation (log-odds) 
%
% the log-odd formula for 2 dimensions is
% x = log(pi/(1-pi))
% for a general diemension,
% x = log(pi/pi(end))
%
% usage:
% x = softmax_transfo_inv(p) operate along the columns 
%
% all the elements  must be in ]0;1[ and the sum along the dimension 
%   must be one.
% 
%
% note: the last row/column in the result equals zero. This make it compatible
% with the reverse function softmax
%
% example:
% --------
% logodds([.9 .1;.2 .8],2)
% softmax(logodds([.9 .1;.2 .8],2),2)
%
%
% See also softmax

x = log(p(1:end-1,:)+1e-320)-repmat(log(p(end,:)+1e-320),size(p,1)-1,1);
