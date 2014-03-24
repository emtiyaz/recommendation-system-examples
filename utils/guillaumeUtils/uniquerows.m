function [urows,v,k] = uniquerows(M)
% UNIQUEROWS - Find unique rows of a matrix 
%
% [urows,ind,ref] = uniquerows(M)
% u:unique rows
% ind:index column vector (length equal to the number of rows)
% ref:reference to one of the unique row (urows = M(ref(i),:))
% See also UNIQUE, matchrows
%
% [a,b,c]=uniquerows(['jj';'kj';'jk';'kk';'kk';'jj'])
[sr,ord] = sortrows(M);
pc = [1;any(diff(sr),2)];
k = ord(find(pc));
urows = M(k,:);
if nargout>1
    cs = cumsum(pc~=0);
    v=zeros(size(M,1),1);
    v(ord) = cs;
end
