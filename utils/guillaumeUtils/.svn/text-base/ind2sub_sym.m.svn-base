function [I,J] = ind2sub_sym(d,idx)
% IND2SUB_SYM - sub-indexes of the unique element in a diagonal matrix
%
%   Example : 
%   [I,J] = ind2sub_sym(3,1:6);[I,J]
% see also find, triu, ones

% [I,J] = find(triu(ones(d)));%slow
[I,J] = cell2lbl(cellfun(@(x) (1:x-1)',num2cell(1:d),'UniformOutput',0));
I = [(1:d)';I];
J = [(1:d)';J];

if nargin>1
    I=I(idx);
    J=J(idx);
end