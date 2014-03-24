function [n,xbar,S2] = first_moments(x,w)

n = sum(w,1);

xbar = x'*w./n;

S2 = (x.^2)'*w;