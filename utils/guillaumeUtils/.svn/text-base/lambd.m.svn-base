function l = lambd(ksi)

%if ksi==0
%    ksi = eps;
%end

I = (abs(ksi)<eps);
ksi(I) = eps;

l = .5./ksi.*(1./(1+exp(-ksi)) - .5);