function I = log_int_exp_quad(A,b,c,inverse)
% log_int_exp_quad - Gaussian integral
% 
% I = log_int_exp_quad(A,b,c)
% compute log of the integral of the following Gaussian function:
%
%   exp(-x'A*x-b'x-c);
%   
% I = log_int_exp_quad(A,b,c,0) %the inverse of A is given
%
% %example: comparison with numerical integration
% 
% a = randn.^2;b=randn;c=randn/3;
% f = @(x) exp(-a*x(:).^2-b*x(:)-c);
% I = log_int_exp_quad(a,b,c)
% x=-10:.01:10; 
% log(trapz(f(x)*diff(x(1:2))))
% plot(x,f(x));
%

if nargin<4 || inverse==0
    A = inv(A);
end

d = length(b);

I = .5*log(det(A)) + .25*b'*A*b - c + d*.5*log(pi);







