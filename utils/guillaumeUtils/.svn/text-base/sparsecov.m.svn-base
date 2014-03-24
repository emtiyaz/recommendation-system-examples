function sigm = sparsecov(X);

nn=size(X,1);
misc.mu = mean(X,1)';
s2bar = X'*X;
sigm = 1./(nn-1)*(s2bar-nn*misc.mu*misc.mu');



