function y = logoneplusexp_bnd(x,xi) 
% LOGONEPLUSEXP_BND - jaakkola's quadratic bound to log(1+exp)
% 
%
% log(1+exp(x)) < L(xi)(x^2-xi^2) + .5*(x-xi) + log(1+exp(xi));
% where L is Jaakkola's function ('lambd' function in statlearn)
% L(xi) = 1/2/xi*(1/(1+exp(x)-.5) 
% 
% %example:
% z=-10:.01:10;plot(z,logoneplusexp(z),'b',z,logoneplusexp_bnd(z,5),'r')
% 
% %bound to the sigmoid
% sigm = @(x) exp(-logoneplusexp(-x));
% sigm_bnd = @(x,xi) exp(-logoneplusexp_bnd(-x,xi));
% z=-10:.01:10;plot(z,sigm(z),'b',z,sigm_bnd(z,5),'r')
%
%
%
% see also logoneplusexp, log1p, lambd, visu_1D_bounds


y = lambd(xi)*(x.^2 - xi.^2) + .5*(x-xi) + logoneplusexp(xi);




