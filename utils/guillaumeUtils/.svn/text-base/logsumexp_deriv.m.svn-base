function [y,d1,d2,d3] = logsumexp_deriv(x)
% Returns log(sum(exp(a),1)) and derivatives
%
% example
% f=@(t)logsumexp_deriv(t);t0=randn(4,1);[D1n,D2n,D3n]=numericderiv(f,t0);[v,D1,D2,D3]=f(t0);[D3n(:) D3(:)],[D2n(:) D2(:)],[D1n(:) D1(:)],[f(t0) v]
%
% see also logoneplusexp

% Written by Guillaume Bouchard


% subtract the largest in each column
[lp,y] = logsoftmax(x,1);
if nargout>1
    d1 = exp(lp);
    if nargout>2
        C = d1*d1';
        d2 = diag(d1) - C;
        if nargout>3
            d = length(x);
            [I,J]= meshgrid(1:d,1:d);
            R = accumarray([I(:) I(:) J(:)],C(:));
            d3 = accumarray([1:d;1:d;1:d]',1:d)...
                + 2*crossmult(d1*d1',d1)...
                - R - permute(R,[3 1 2]) - permute(R,[2 3 1]);
            if d==2%correct some strange bug for d=2...
                d3(1,1,1) = -d3(2,1,1);
                d3(2,2,2) = -d3(1,2,2);
            end                
        end
    end
end