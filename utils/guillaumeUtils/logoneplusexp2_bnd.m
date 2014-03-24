function y = logoneplusexp2_bnd(a,x,xi) 
% LOGONEPLUSEXP_BND - Jebara's quadratic bound to log(a+exp(x^2))
% 
%
% log(a+exp(x^2/2)) < log(a+exp(xi^2/2)) + 1/(1+a*exp(-xi^2/2)xi(x-xi)
%                     + 1/2*(x-xi)^2*(1+J(a,xi)*xi^2);
% where J is Jebara's function ('jebara' function in statlearn)
% J(xi) = tanh(.5*(log(a)-.5*xi^2))/2/(log(a)-.5*xi^2) 
%  
% %example:
% a=randn^2;xi=randn*3;x=linspace(-25,25,1000);f=@(x)log(a+exp(x.^2/2));plot(x,f(x),'b',x,logoneplusexp2_bnd(a,x,xi),'r',[xi xi],[0 f(xi)],'k.-');
%
% see also logoneplusexp, log1p, lambd, visu_1D_bounds

xi2 = xi.^2;
la=log(a);
temp = .5*xi2-la;
y0 = la + logoneplusexp(temp);
dy0 = exp(-logoneplusexp(-temp))*xi;

y = y0 + (x-xi)*dy0 + (x-xi).^2/2*(1+jebara(a,xi)*xi2); 



