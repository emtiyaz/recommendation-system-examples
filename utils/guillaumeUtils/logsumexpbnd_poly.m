function [a,b,c] = logsumexpbnd_poly(lvp,mat)
% logsumexpbnd_poly - coefficient of the quadratic upper bound to the log-sum-exp funciton
% [a,b,c] = logsumexpbnd_poly(lvp)
% the bound is x^2*a + b*x + c
% in matlab notation, the bound is 
% sum(x.^2.*(ones(size(x,1),1)*a),2) + x*b' + c
% 
% [a,b,c] = logsumexpbnd_poly(lvp,0) return horizontal vectors 'a' and 'b' (default)
% [a,b,c] = logsumexpbnd_poly(lvp,1) return a matrix for 'a'
%
% see logsumexpbnd_optim, logsumexpbnd_optim_demo

if nargin<2
    mat=0;
end
m = lvp.method;
if ischar(lvp.method) && strcmp('worseCurvature',lvp.method)
    m=0;
end




switch m
    case 'taylorOnMean' %laplace (matrix)
        [n,K] = size(lvp.xi);
        K = length(lvp.xi);        
        [s,ls] = softmax(lvp.xi,2);
        a = .5*(diag(s) - s'*s);
        b = (s' - 2*a*lvp.xi')';        
        c = lvp.xi*a*lvp.xi' - lvp.xi*s' + ls;
    case -1%diagonal bound (matrix)
        [n,K] = size(lvp.xi);
        K = length(lvp.xi);
        a = eye(K)/K/2;
        [s,ls] = softmax(lvp.xi,2);
        b = (s' - 2*a*lvp.xi')';
        c = lvp.xi*a*lvp.xi' - lvp.xi*s' + ls;
    case 0%bohnings bound (matrix)
        [n,K] = size(lvp.touchPoint);
        K = length(lvp.touchPoint);
        a = (eye(K)/K - ones(K)/K^2)/2;
        %a = eye(K)/2;
        [s,ls] = softmax(lvp.touchPoint,2);
        b = (s' - 2*a*lvp.touchPoint')';
        c = lvp.touchPoint*a*lvp.touchPoint' - lvp.touchPoint*s' + ls;
    case 1
        [n,K] = size(lvp.xi);
        a = lambd(lvp.xi);
        b = .5-2*(lvp.scale*ones(1,K)).*a;
        c = lvp.scale.*(1-K/2 + sum(a,2).*lvp.scale) - sum(lvp.xi.*(a.*lvp.xi+.5),2) + sum(logoneplusexp(lvp.xi),2);
        
    case 2
        [n,K] = size(lvp.xi);
        [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi',0);
        a = lambd(lvp.xi)./lvp.epsi;        
        b = 1./lvp.epsi*(.5-(1-lvp.epsi)*df)';
        cu = sum(logoneplusexp(lvp.xi)-lambd(lvp.xi).*lvp.xi.^2-lvp.xi/2);
        cv = f-lvp.chi*df;      
        c = (cu - (1-lvp.epsi)*cv - entrop(lvp.epsi))./lvp.epsi;
        
        
    case {'simpleLooseBound','tightBoundNoScale','tightBound'}
        [n,K] = size(lvp.xi);
        [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi(:),0);
        lb = lambd(lvp.xi(:))';
        a = lb./lvp.epsi;
        b = 1./lvp.epsi*(.5-2*(lvp.scale*ones(1,K)).*lb-(1-lvp.epsi)*df');
        cu = sum(logoneplusexp(lvp.xi)+lb.*(lvp.scale.^2-lvp.xi.^2) - (lvp.xi + lvp.scale)/2);
        cv = f-(lvp.chi+lvp.scale)*df;
        c = lvp.scale + (cu - (1-lvp.epsi)*cv - entrop(lvp.epsi))./(lvp.epsi+(lvp.epsi==0));
        
    case {'tightBound2'}
        [n,K] = size(lvp.xi);
        [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi(:),0);
        lb = lambd(lvp.xi(:))';
        a = lb./lvp.epsi;
        b = 1./lvp.epsi*(.5-2*(lvp.scale*ones(1,K)).*lb-(1-lvp.epsi)*df');
        cu = sum(logoneplusexp(lvp.xi)+lb.*(lvp.scale.^2-lvp.xi.^2) - (lvp.xi + lvp.scale)/2);
        cv = f-(lvp.chi+lvp.scale)*df;
        c = lvp.scale + (cu - (1-lvp.epsi)*cv - entrop(lvp.epsi))./(lvp.epsi+(lvp.epsi==0));

        
        
        
    case 4  
        [n,K] = size(lvp.xi);
        if length(lvp.epsi)==1
            lvp.epsi(2)=1-lvp.epsi(1);
        end
        [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi(:),1);
        lb = lambd(lvp.xi(:))';
        a = lb./lvp.epsi(1);
        b = 1./lvp.epsi(1)*(.5-2*(lvp.scale*ones(1,K)).*lb-lvp.epsi(2)*df');
        cu = sum(logoneplusexp(lvp.xi)+lb.*(lvp.scale.^2-lvp.xi.^2) - (lvp.xi + lvp.scale)/2);
        cv = f-(lvp.chi+lvp.scale)*df;
        c = lvp.scale + (cu - lvp.epsi(2)*cv - myentrop([lvp.epsi 1-sum(lvp.epsi)]'))./(lvp.epsi(1)+(lvp.epsi(1)==0));
    case 'minTightBound'
        if lvp.choice==0
            lvp.method = 'worseCurvature';
            [a,b,c] = logsumexpbnd_poly(lvp,mat);
        else
            lvp.method = 'tightBound';            
            [a,b,c] = logsumexpbnd_poly(lvp,mat);
        end
end

if mat 
    if min(size(a))==1%not matrix
        a = diag(a);        
    end
else%vector
    if min(size(a))>1%matrix
        a = diag(a)';
    end
end

if isfield(lvp,'beta') %bound improved by a linear transfom    
    T = eye(K)-ones(K,1)*lvp.beta';
    a = T'*a*T;%assumes a is a matrix
    b = (T'*b' + lvp.beta)';    
end


if nargout==1
    a = {a,b,c};
end

function e=myentrop(x)
e=sum(x.*log(x+(x==0)));
