function [I] = matchrows(A,B)
% MATCHROWS - match rows of 2 matrices
%
% I = matchrows(A,B)
%   returns a vector I of length size(A,1) such that A(i,:) ==  B(I(i),:)
%   or 0 if there is no row in B equal to A(i,:)
%   B should contain only unique rows
%
% example:
% A=['jj';'kj';'jk';'kk';'kk';'jj'];
% B=['kj';'kk';'jj'];
% I = matchrows(A,B)
%

nA = size(A,1);
nB = size(B,1);
binvec = [zeros(nA,1);ones(nB,1)];
[U,ctg,sel] = uniquerows([A;B]);
ctgA = ctg(1:nA);
ctgB = ctg(nA+1:end);
map = zeros(1,max(ctgA));
map(ctgB) = 1:nB;
I =  map(ctg(1:nA));

% [sr,ord] = sortrows([A;B]);
% Ichg = 
% cumsum([1;diff(binvec(ord))>0])
% ord(Ichg)

% pc = [1;any(diff(sr),2)];
% k = ord(find(pc));
% urows = M(k,:);
% if nargout>1
%     cs = cumsum(pc~=0);
%     v=zeros(size(M,1),1);
%     v(ord) = cs;
% end
