function y = logoneplusexpabs_bnd(x,xi0,a,b) 
% LOGONEPLUSEXPABS_BND - quadratic bound to log(1+exp(a*abs(x)+b)
% 
% L is Jaakkola's function: L(xi) = 1/2/xi*(1/(1+exp(x)-.5) 
% 
% if L(xi)>1/4b
%    y = Lxi*(x.^2 + b^2 - xi.^2) + (.5-2*b*Lxi)*x - .5*xi0 + logoneplusexp(xi);
% otherwise,
%     y = Lxi*(x.^2 + b^2 - xi.^2) + (.5-2*b*Lxi)*((x.^2)/2./xi0 + xi0/2) - .5*xi0 + logoneplusexp(xi);
%
% %examples:
% %--------
% figure(1)
% a=-1;b=3;f = {@(x) log(1+exp(a*abs(x)+b))};for xi=[-3 1 8]; f{end+1} = @(x)(logoneplusexpabs_bnd(x,xi,a,b));end;visu_1D_bounds(f,'lims',[-10 10]);
% figure(2)
% a=1;b=3;f = {@(x) log(1+exp(a*abs(x)+b))};for xi=[-3 1 8]; f{end+1} = @(x)(logoneplusexpabs_bnd(x,xi,a,b));end;visu_1D_bounds(f,'lims',[-10 10]);
% figure(3)
% a=1;b=-3;f = {@(x) log(1+exp(a*abs(x)+b))};for xi=[-3 1 8]; f{end+1} = @(x)(logoneplusexpabs_bnd(x,xi,a,b));end;visu_1D_bounds(f,'lims',[-10 10]);
% 
% see also logoneplusexp, logoneplusexp_bnd, log1p, lambd, visu_1D_bounds

xi = a*abs(xi0) + b;


Lxi = lambd(xi);
x=a*x;    

if a<0 || 4*Lxi*b<-1      
    y = Lxi*(x.^2 + b^2 - xi.^2) + sign(xi0)*(.5+2*b*Lxi)*x - .5*(xi-b) + logoneplusexp(xi);
else
    y = Lxi*(x.^2 + b^2 - xi.^2) + (.5+2*b*Lxi)*((x.^2)/2./abs(xi-b) + abs(xi-b)/2) - .5*(xi-b) + logoneplusexp(xi);            
end    



