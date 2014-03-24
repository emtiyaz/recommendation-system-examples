function decs = pack_state(bins)
% PACK_STATE - transform the full state space variable into a unique value 
% 
% x = pack_state(S,x);
%
% example:
% --------
%
% pack_state(unpack_state(1:8,3))
%
% See also unpack_state


% bins = bins(:,end:-1:1);

decs = bin2dec(char(bins(:,1:end)+'0'-1))+1;


