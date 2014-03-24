function y = log_exp_min_one(x) 
% log_exp_min_one
%
% log((exp(x)-1)) may have some overflow. This function prevents it.
%
% x=logspace(-4,4,1000);semilogx(x,log((exp(x)-1)),'b.',x,log_exp_min_one(x),'r')

big=x>10;

y = zeros(size(x));
if any(~big)
    %y(~big) = log(exp(x(~big))-1);
%     %MODIFIED BY EMTIYAZ
%     % exp(x) can become less than 1, which gives imaginary values, I reset
%     % these value to be 1+eps
%     t = exp(x(~big));
%     idx = find(t<1);
%     t(idx) = 1 + eps;
%     ynew = log(t-1);
    y(~big) = log(max(1e-300,exp(x(~big))-1));
end
if any(big)
    y(big) = x(big) + log(1-exp(-x(big)));
end

