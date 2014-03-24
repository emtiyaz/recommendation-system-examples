function bins = unpack_state(decs,ndims)
% UNPACK_STATE - transform factored space into full state space variable
% 
% S = unpack_state(c,x);
%
% example:
% --------
%
% pack_state(unpack_state(1:8,3))
%
% See also pack_state
% if nargin>1
    bins = (dec2bin(decs-1)-'0'+1);
    bins = [ones(length(decs),ndims -size(bins,2)) bins];
% else
%     bins = (dec2bin(decs-1)-'0'+1);
% end
bins = bins(:,1:end);
