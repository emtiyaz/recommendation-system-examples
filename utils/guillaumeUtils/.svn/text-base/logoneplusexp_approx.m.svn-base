function [y,a,b,c,gm] = logoneplusexp_approx(x,m,taylor2) 
% LOGONEPLUSEXP_APPROX - Gaussian + max approximation to log(1+exp)
% 
% logoneplusexp_approx(x) 
% makes a Taylor expansion of log(log(1+exp(x))) at x=0 valid for negative x
%
% logoneplusexp_approx(x,m)
% makes a 1st order Taylor expansion of log(log(1+exp(x))) around x=-abs(m)
% and symetrize
%
% [y,a,b,c] = logoneplusexp_approx(x,m)  also return coefficients so that
% the function equals y = exp(-.5*a*x.^2+b*abs(x)+c) + max(x,0);
%
% %example:
% z=-10:.01:10;plot(z,logoneplusexp(z)-logoneplusexp_approx(z),'r')
% m=randn*10;z=-10:.01:10;plot(z,logoneplusexp(z)-logoneplusexp_approx(z,m),'r',m,0,'.')
% m=-.2;z=m-2:.001:m+2;plot(z(2:end-1),diff([logoneplusexp(z);logoneplusexp_approx(z,m)],2,2));
% 
% 
% %example of derivative according to the approximation point
% z=rand*5;m=-10:.01:10;[f,a,b,c,gm]=logoneplusexp_approx(z,m);plot(m,gm,'b',m(2:end),diff(f)./diff(m(1:2)),'r--')
% 
% see also logoneplusexp, log1p, lambd, visu_1D_bounds


if nargin<2
    
    L2=log(2);
    a = (1-L2)/4/L2;
    b = -1/2/L2;
    c = log(L2);
        
    y = exp(-.5*a*x.^2+b*abs(x)+c) + max(x,0);
else
    si = sign(m);
    m = -abs(m);
    L = logoneplusexp(m)+1e-300;
    lL = log(L);
    sigmm = exp(m-L);
    sigmmc = exp(-L);
    d1 = sigmm./L;
    if nargin<3 || taylor2
        d2 = d1.*(sigmmc-d1);        
    else
        d2 = 2*(log(log(2))-lL+d1.*m)./m.^2;
    end
%     y0 = L*exp(d1*(-abs(x)-m)+.5*d2*(-abs(x)-m).^2)+ max(x,0);
    
    a = -d2;
    b = -d1+m.*d2;
    c = lL - d1.*m + d2.*m.^2/2;
    
    expe = exp(-.5*a.*x.^2+b.*abs(x)+c);
    y = expe + max(x,0);

end

if nargout>4
    dL = sigmm;
    dsigmm = sigmmc.*sigmm;
    dsigmmc = -dsigmm;
    dd1 = (dsigmm.*L-sigmm.*dL)./L.^2;
    dd2 = dd1.*(sigmmc-d1) + d1.*(dsigmmc-dd1);
    da = -dd2;
    db = -dd1 + d2 + m.*dd2;
    dc = dL./L - dd1.*m - d1 + dd2.*m.^2/2 + d2.*m;
    gm = (-.5*da.*x.^2+db.*abs(x)+dc).*expe;
    
    gm = -gm.*(si+(m==0));%recover the signs
end
