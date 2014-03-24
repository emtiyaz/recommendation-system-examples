function lvp = logsumexpbnd_optim(mu,var,lvp,maxiter,method)
% logsumexpbnd_optim - optimize the log-sum-exp bound
%
% usage:
% lvp = logsumexpbnd_optim(mu,var)
% lvp = logsumexpbnd_optim(mu,var,old_lvp)
% lvp = logsumexpbnd_optim(mu,var,old_lvp,maxiter)
% lvp = logsumexpbnd_optim(mu,var,old_lvp,maxiter,method)
%
% Example:
%
% beta=[-10 -20;0 0;3 0;10 -20];
% fplot(@(x) logsumexp(x*beta(:,1)+beta(:,2),1),[-5,5]);
% mu = 1;var=.001;
% lvp=logsumexpbnd_optim(mu*beta(:,1)'+beta(:,2)',beta(:,1)'.^2*var);
% [a,b,c]=logsumexpbnd_poly(lvp);hold on;
% fplot(@(x) sum((x*beta(:,1)'+beta(:,2)').^2.*a,2) + sum((x*beta(:,1)'+beta(:,2)').*b,2) + c,[-5,5],'r');
% set(gca,'YLim',[0 30]);

%% demo

if nargin<1%demo
    %calls a script
    logsumexpbnd_optim_demo
end

%% main program
[n,K] = size(mu);
if  nargin<3
    lvp = [];
end

if nargin<4 || isempty(maxiter);
    maxiter = 100;
end

if nargin<5
    method = 'tightBound';
end

if min(size(var))>1
    var = diag(var)';
end

switch method
    %% Bohning's bound
    case 'taylorOnMean'
        lvp.xi = mu;       
    case 'Bohning'
        lvp.xi = mu;       
        %% Bohning's bound
    case 'worseCurvature'
        lvp.touchPoint = mu;       
        
    case 'simpleLooseBound'
        lvp = lvp_update(mu,var,lvp,[1,1,0]);
        
    case 'tightBoundNoScale'
        lvp = lvp_update(mu,var,lvp,[1,0,1]);
        
    case 'tightBound'
        lvp = lvp_update(mu,var,lvp,[1,1,1]);
    case 'tightBound2'
        lvp = lvp_update3(mu,var,lvp,[1,1,1]);
        
    case 'minTightBound'
        lvp = lvp_update(mu,var,lvp,[1,1,1]);
        lvp.touchPoint = mu;       

        lvp.method = 'worseCurvature';
        [A0,b0,c0] = logsumexpbnd_poly(lvp,1);
        lvp.method = 'tightBound';            
        [A1,b1,c1] = logsumexpbnd_poly(lvp,1);
        
        lvp.choice = mu*A1*mu' + mu*b1' + c1 + var*diag(A1)...
            < mu*A0*mu' + mu*b0' + c0 + var*diag(A0);
        
%     case 'tightBound2'
%         lvp = lvp_update2(mu,var,lvp,[1,1,1,1]);
       
    otherwise
        error('invalid method')
end
lvp.method = method;



%% utility function solving a simple problem involving entropy
function x = upd_min_y_minus_Hx_ov_x(y,x0)
% 
% epsil = 1e-5;
% func = @(x) (y + x*log(max(x,epsil)) + (1-x)*log(max(1-x,epsil)))/max(x,epsil);
% 
% opts = optimset('TolX',1e-12,'Display','off','MaxIter',100);
% if exp(y) < 1-x0
%     x = fminbnd(func,0,x0,opts);
% else
%     maxx=1-epsil;
%     x = fminbnd(func,min(maxx,x0),maxx,opts);
% end
% 
% if func(x0)<func(x) %set old value if no improvement
%     x = x0;
% end

x = 1-exp(-y);

%OLD STUFF BASED ON THE QUADRATIC LOWER BOUND TO THE ENTROPY THAT CONVERGES TO SLOWLY TO 0 OR 1
%to view the algorithm iterations, launch:
%x=linspace(.01,1-eps,10000);clf;D=.1;z=.7;for i=1:6;[a,b,c]=entropy_bnd(x,z);line(x,1./x.*(D-a*x.^2-b*x-c),'Color',[i/6 1 0]);z=sqrt(-(D-c)/a);line(z,1./z.*(D-a*z.^2-b*z-c),'Marker','+','Color','k');end;line(x,1./x.*(D-entrop(x)),'LineWidth',2);set(gca,'YLim',[-3 2]);xlabel('\epsilon');ylabel('\epsilon^{-1}(D-H(\epsilon))');title(['D=' num2str(D)])
%[a,b,c]=entropy_bnd(x0);
%x = min(sqrt(-(y-c)/a),1-eps);

function [x,y] = upd_u_min_yv_minus_Hxy_ov_x(u,v,x0,y0)
epsil = 1e-15;
opts = optimset('TolX',1e-12,'Display','off','MaxIter',10);
%% update x (same as before)
C = u-y0*v + y0*log(max(y0,epsil));
func = @(x) (C + x*log(max(x,epsil)) + (1-x-y0)*log(max(1-x-y0,epsil)))/max(x,epsil);
if exp(C) < 1-x0
    x = fminbnd(func,0,x0,opts);
else
    maxx=1-y;
    x = fminbnd(func,min(maxx,x0),maxx,opts);
end
if func(x0)<func(x) %set old value if no improvement
    x = x0;
end

%% update y
funcy = @(y) u - y*v + x*log(max(x,epsil)) + y*log(max(y,epsil) + (1-x-y)*log(max(1-x-y,epsil)));
[y1,f1] = fminbnd(funcy,0,y0,opts);%left
maxy = 1 - x;
[y2,f2] = fminbnd(funcy,min(maxy,y0),maxx,opts);%right
if f1<f2
    y=y1;
else
    y=y2;
end
if funcy(y0)<funcy(y) %set old value if no improvement
    y = y0;
end


function R = xMx(x,M)
R = x*M*x';

%#ok<*DEFNU>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,fnew,rate] = upd_min_entropy_subfunc(x0,y0,a,b,rate)
eps = 1e-5;
eps2 = 1e-10;
x0 = max(eps,x0);
y0 = max(eps,y0);
if x0+y0>=1
    x0 = (.5-eps)*(x0-y0);
    y0 = (.5-eps)*(y0-x0);
end
[f,df] = entropy_subfunc(x0,y0,a,b);
x=x0;
y=y0;
fnew = inf;
while fnew>f
    rate = rate*.9;
    df = df./((df==0)+norm(df));
    x = max(eps,min(1-eps2,x0-df(1)*rate));
    y = max(eps,min(1-eps2,y0-df(2)*rate));
    s = x+y;
    if s>1-eps2
        x = max(eps2,x-s/2);
        y = max(eps2,y-s/2);
    end
    [fnew,df] = entropy_subfunc(x,y,a,b);
end

rate = rate*1.5 + .1e-3;



%% check validity of updates
function check_update(lvp_old,lvp,mu,var,msg)

if ~isempty(lvp_old) && isfield(lvp_old,'chi')
    [a_old,b_old,c_old] = logsumexpbnd_poly(lvp_old,1);
    E_old = xMx(mu,a_old) + var*diag(a_old) + mu*b_old' + c_old;
else
    E_old = inf;
end
if ~isempty(lvp) && isfield(lvp,'chi')
    [a,b,c] = logsumexpbnd_poly(lvp,1);
    E = xMx(mu,a) + var*diag(a) + mu*b' + c;
    
    fprintf('%20s: %f (%5f)\n',msg,E,E_old-E);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lvp = lvp_update(mu,var,lvp,scheme)

K = length(mu);
if isempty(lvp)
    if scheme(3)
        lvp.epsi = .5;
    else
        lvp.epsi = 1;
        lvp.epsi0 = 1;
    end
    if scheme(2)
        lvp.scale = max(mu);
    else
        lvp.scale = 0;
    end
    if scheme(1)
        lvp.xi = mu;
        lvp.chi = mu;
    else
        lvp.xi = zeros(1,K);
        lvp.chi = zeros(1,K);
    end
end

if scheme(1) %xi updates
    E = mu - lvp.scale;
    E2 = E.^2 + var;
    lvp.xi = sqrt(E2);
end


if scheme(3) % epsi updates
    lvp.chi = mu - lvp.scale;
    
    if 1
      lb = lambd(lvp.xi);
      u = sum( lb.*(E2-lvp.xi.^2) + .5*(E - lvp.xi) + logoneplusexp(lvp.xi) );
      [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi',0);
      v = f + (E-lvp.chi)*df;  
  %   lvp_old = lvp;
      lvp.epsi = upd_min_y_minus_Hx_ov_x(u-v,lvp.epsi);

    else
      %Written by Emt
      lvpTemp = lvp;
      lvpTemp.method = 'simpleLooseBound';
      lvpTemp.epsi = 1;
      [A,b,c] = logsumexpbnd_poly(lvpTemp,1);
      f1 = trace(A*(mu(:)*mu(:)'+diag(var(:)))) + b(:)'*mu(:) + c;
      [C_chi,g_chi] = log_prod_1_plus_exp_min_sum_exp(lvp.chi',0);
      f2 = C_chi - g_chi(:)'*(lvp.chi(:) + lvp.scale -mu(:)) + lvp.scale;
      lvp.epsi = 1-exp(f2-f1);
    end
      
%     %check convergence
%     lvp_old.method = 3;
%     lvp.method = 3;
%     check_update(lvp_old,lvp,mu,var,'epsi')
%     lvp_old = lvp;
%     
end



if scheme(2) %scale updates
    lb = lambd(lvp.xi);
    [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi',0);
    
    lvp.scale = (sum(mu.*lb - (1-lvp.epsi)/2*df',2) + (K/2-lvp.epsi)/2)./(sum(lb,2)+eps);
    s=lvp.scale;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lvp = lvp_update2(mu,var,lvp,scheme)
K = length(mu);
if isempty(lvp)
    if scheme(3)
        lvp.epsi = [.5 .5];
    else
        lvp.epsi = [1 0];
    end
    if scheme(4)
        lvp.epsi = [.5 .4];
    end
    if scheme(2)
        lvp.scale = max(mu);
    else
        lvp.scale = 0;
    end
    if scheme(1)
        lvp.xi = mu;
        lvp.chi = mu;
    else
        lvp.xi = zeros(1,K);
        lvp.chi = zeros(1,K);
    end
end

if scheme(1) %xi updates
    E = mu - lvp.scale;
    E2 = E.^2 + var;
    lvp.xi = sqrt(E2);
end


if scheme(3) % epsi updates
    lvp.chi = mu - lvp.scale;
    
    lb = lambd(lvp.xi);
    u = sum( lb.*(E2-lvp.xi.^2) + .5*(E - lvp.xi) + logoneplusexp(lvp.xi) );
    [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi',1);
    v = f + (E-lvp.chi)*df;  


    
    lvp.epsi(1) = upd_min_y_minus_Hx_ov_x(u-v,lvp.epsi(1));
    lvp.epsi(2) = 1-lvp.epsi(1);
    
%     
end

%      lvp_old = lvp;
if scheme(2) %scale updates
    lb = lambd(lvp.xi);
    [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi',1);
    lvp.scale = (sum(mu.*lb - lvp.epsi(2)/2*df',2) + (K/2-lvp.epsi(1))/2)./(sum(lb,2)+eps);
end

%     %check convergence
%     lvp_old.method = 4;
%     lvp.method = 4;
%     check_update(lvp_old,lvp,mu,var,'epsi')
%     lvp_old = lvp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lvp = lvp_update3(mu,var,lvp,scheme)

K = length(mu);
if isempty(lvp)
    if scheme(3)
        lvp.epsi = .5;
    else
        lvp.epsi = 1;
        lvp.epsi0 = 1;
    end
    if scheme(2)
        lvp.scale = max(mu);
    else
        lvp.scale = 0;
    end
    if scheme(1)
        lvp.xi = mu;
        lvp.chi = mu;
    else
        lvp.xi = zeros(1,K);
        lvp.chi = zeros(1,K);
    end
end

if scheme(1) %xi updates
    E = mu - lvp.scale;
    E2 = E.^2 + var;
    lvp.xi = sqrt(E2);
end


if scheme(3) % epsi updates
    lvp.chi = mu - lvp.scale;
    
    lb = lambd(lvp.xi);
    u = sum( lb.*(E2-lvp.xi.^2) + .5*(E - lvp.xi) + logoneplusexp(lvp.xi) );
    [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi',1);
    v = f + (E-lvp.chi)*df;  

%     lvp_old = lvp;
    
    lvp.epsi = upd_min_y_minus_Hx_ov_x(u-v,lvp.epsi);
    
%     %check convergence
%     lvp_old.method = 3;
%     lvp.method = 3;
%     check_update(lvp_old,lvp,mu,var,'epsi')
%     lvp_old = lvp;
%     
end



if scheme(2) %scale updates
    lb = lambd(lvp.xi);
    [f,df] = log_prod_1_plus_exp_min_sum_exp(lvp.chi',1);
    
    lvp.scale = (sum(mu.*lb - (1-lvp.epsi)/2*df',2) + (K/2-lvp.epsi)/2)./(sum(lb,2)+eps);
    s=lvp.scale;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
