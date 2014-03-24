function y = logexpminusoneovx(x) 
% LOGEXPMINUSONEOVX - Secure logarithm of (exp(x)-1)/x
%
% log((exp(x)-1)/x) may have some overflow. This function prevents it.
%
% useful to compute the log of the 'lambd' function used to bound the sigmoid:
% %example (compare with 'log(lambd(x))')
% x = linspace(-20,20,1000);plot(x,log(lambd(x)),'b',x,log(.25)+logexpminusoneovx(x)-logoneplusexp(x),'r--')
%
% Remark: the first derivative is exp(x-f(x))/x
% %example(derivative)
% x = linspace(-20,20,1000);
% plot(x,(exp(x-logexpminusoneovx(x))-1)./x,'b',x(2:end),diff(logexpminusoneovx(x))/diff(x(1:2)),'r')
%
% %example (demonstrate the overflow):
% log(exp(1000)-1)/1000
% logexpminusoneovx(1000)
% 
% example 2: the function plot
% x=linspace(-100,100,10000);plot(x,log((exp(x)-1)./x),'b',x,logexpminusoneovx(x),'r')
%
% %example 2 (plot the function at transition points)
% x = linspace(-.01,.01,1000)-33;plot(x(2:end-1),diff(logexpminusoneovx(x),2));
% x = linspace(-.01,.01,1000);plot(x(2:end-1),diff(logexpminusoneovx(x),2));
%
%
% see also logsumexp

if nargin<2
    dim=1;
end
thr1 = -33;
thr2 = 33;

I = 1+(x>thr1)+(x>=0)+(x>0)+(x>thr2);
C = lbl2cell(I,[],'K',5);
y = zeros(size(x));

y(C{1}) = - log(-x(C{1}));%large negative
% y(C{2}) = log(1-exp(x(C{2}))) - log(-x(C{2}));%negative
y(C{2}) = log((exp(x(C{2}))-1)./x(C{2}));
y(C{3}) = 0;
% y(C{4}) = log(exp(x(C{4}))-1) - log(x(C{4}));
y(C{4}) = log((exp(x(C{4}))-1)./x(C{4}));
y(C{5}) = x(C{5})-log(x(C{5}));%asymptot large positive



