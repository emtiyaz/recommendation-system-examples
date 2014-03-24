function x = logodds(p,dim) 
% LOGODDS - Log-odds function on proportions (inverse softmax)
%
% the log-odd formula for 2 dimensions is
% x = log(pi/(1-pi))
% for a general diemension,
% x = log(pi/pi(end))
%
% usage:
% x = logodds(p) %along the first dimension (columns)
% x = logodds(p,2) %along the rows
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

if nargin<2
    dim=1;
end
if dim==1
    x = log(p(1:end-1,:)+1e-320)-repmat(log(p(end,:)+1e-320),size(p,1)-1,1);
elseif dim==2
    x = log(p(:,1:end-1)+1e-320)-repmat(log(p(:,end)+1e-320),1,size(p,1)-1);
end
