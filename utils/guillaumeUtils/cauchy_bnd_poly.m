function [a,b,c] = cauchy_bnd_poly(xi)
% cauchy_bnd_poly - coefficient of the quadratic lower bound to the (log-)Cauchy pdf
% [a,b,c] = cauchy_bnd_poly(lvp)
% the bound is 1./(1+x.^2) >= -.5*x^2*a(xi) + b(xi)*x + c(xi) 
% in matlab notation, the bound is 
%   -.5*x.^2*a + x*b + c
% 
% %example: 
% x=-10:.01:10;
% xi=randn;[a,b,c]=cauchy_bnd_poly(xi);plot(x,1./pi./(1+x.^2),'b',x,exp(-.5*x.^2*a + x*b + c),'r')
%
% see logsumexpbnd_poly


a = 2./(1+xi^2);
b = 0;
c = xi^2/(1+xi^2)-log(pi*(1+xi^2));

if nargout<=1
    a = [a;b;c];
end


