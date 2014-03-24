function [y,grad,d2,d3] = logoneplusexp(x) 
% LOGONEPLUSEXP - Secure logarithm of (1+exp(x))
%
% log(1+exp(x)) may have some overflow. This function prevents it.
% 
% %example (demonstrate the overflow):
% log(1+exp(1000))
% logoneplusexp(1000)
%
% %example 2 (plot the function at transition points)
% x = linspace(-17-.01,-17+.01,1000);plot(x(2:end-1),diff(logoneplusexp(x),2));
% x = linspace(-.01,.01,1000);plot(x(2:end-1),diff(logoneplusexp(x),2));
% x = linspace(33-.01,33+.01,1000);plot(x(2:end-1),diff(logoneplusexp(x),2));
% 
% % demo of derivatives:
% x=linspace(-5,5,10000);[f,g]=logoneplusexp(x);plot(x(2:2:end-1),[diff(f(1:2:end))/2/diff(x(1:2))-g(2:2:end-1)]')
% 
% % approximation using Gaussian + max
% plot(x,[logoneplusexp(x)-(exp(-.5*a*x.^2+b*abs(x)+c)+max(x,0))])
%
% %third derivative
%  f=@(t) logoneplusexp(t);t0=randn(1,1);[D1n,D2n,D3n]=numericderiv(f,t0);[v,D1,D2,D3]=f(t0);[D3n(:) D3(:)],[D2n(:) D2(:)],[D1n(:) D1(:)]
% see also logsumexp

% % max(x,0) could be used as an approximation:
% y = max(0,x);
% return


y0 = log(1+exp(-abs(x)));
y = y0 + max(x,0);

   
if nargout>1
    s = exp(-y);
    grad = 1-s;
    if nargout>2      
        d2 = s.*grad;
        if nargout>2
            d3 = d2.*(2*s-1);
        end
    end
end



if 0
    
    thr1 = -17;
    thr2 = 33;
    % cell array
    % I = 1+(x>thr1)+(x>0)+(x>thr2);
    % C = lbl2cell(I,[],'K',4);
    
    C = findintervals(x,[thr1 0 thr2],1);
    
    y = zeros(size(x));
    y(C{1}) = exp(x(C{1}));%asymptot large negative
    ex = exp(x(C{2}));
    l1p = 1+ex;
    y(C{2}) = log(l1p);%negative
    exm = exp(-x(C{3}));
    l1pm = 1+exm;
    y(C{3}) = x(C{3}) + log(l1pm);%positive
    y(C{4}) = x(C{4});%asymptot large positive
    
    if nargout>1
        %gradient:sigmoid
        grad = zeros(size(x));
        
        grad(C{1}) = exp(x(C{1}));%asymptot large negative
        grad(C{2}) = ex./l1p;%negative
        grad(C{3}) = 1./l1pm;%positive
        grad(C{4}) = 1-exp(-x(C{4}));%asymptot large positive
        
    end
end


