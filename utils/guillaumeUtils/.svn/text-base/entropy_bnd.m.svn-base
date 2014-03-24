function [a,b,c] = entropy_bnd(x0,x)
% entropy_bnd - quadratic lower bound to the entropy
% 
% % examples
%   entropy_bnd     %demo
%   entropy_bnd(x0) %demo for approximatin at x0
%
% %plot example
%   x0=rand;y1=-x.*log(x) - (1-x).*log(1-x);y2 = entropy_bnd(x0,x);plot(x,y1,'b',x,y2,'r',x0,entropy_bnd(x0,x0),'k.');
%
% see also entrop, logsumexpbnd_optim

if nargin<2 && nargout<=1
    x = linspace(.001,.999,1000);
    if nargin<1
        x0 = rand;
    end
    y1 = -x.*log(x) - (1-x).*log(1-x);
    y2 = entropy_bnd(x0,x);
    plot(x,y1,'b',x,y2,'r',x0,entrop(x0),'k.');
    grid
    set(gca,'YLim',[-1,1])
    return
end

l1 = log(x0);
l2 = log(1-x0);
a = max(l1,l2)/min(x0,1-x0).^2;
c = - l2 + x0.^2*a;
b = -l1 + l2 - 2*x0*a;

if nargout<=1
    a = a*x.^2 + b*x + c;
end
