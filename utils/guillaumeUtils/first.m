function Y=first(X,dim)
% FIRST - first non-zeros elements
%
% Y=first(X) return the first non-zero element of the vector X
% 
% For matrices,
%   Y=first(X,1) [default] operates along the columns
%   Y=first(X,2) operates along the rows
%
% Note: return 0 if the whole row/column is null
%
% example:
% --------
%
% first([0 0 -10 3 0;0 19 12 0 0],2)
% first([0 0 -10 3 0;0 19 12 0 0],1)



if nargin>=2 & (dim==2)
    X=X';
end

[n,d] = size(X);
[null,I] = max(X~=0);
Y = X(sub2ind([n,d],I,1:d));

