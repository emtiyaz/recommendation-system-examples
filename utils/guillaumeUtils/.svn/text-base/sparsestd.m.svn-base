function sigm = sparsestd(X);
% sparsestd - standard deviation for sparse matrix
%
% see also sparsecov

nn = size(X,1);
mu = (mean(X,1));
s2bar = (sum(X.^2,1));
sigm = sqrt(1./(nn-1)*(s2bar-nn*mu.^2));




