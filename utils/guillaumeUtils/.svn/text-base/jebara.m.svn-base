function y = jebara(a,x)
% jebara - curvature of the bound to the logistic-quadratic function    
% utility function used in logoneplusexp2_bnd
%
% % example:
% a=.01;x=linspace(-15,15,1000);y=jebara(a,x);plot(x,y);
x2 = x.^2;

lamhx2 = log(a)-.5*x2;

y = tanh(.5*lamhx2)./2./(lamhx2+(lamhx2==0)); 