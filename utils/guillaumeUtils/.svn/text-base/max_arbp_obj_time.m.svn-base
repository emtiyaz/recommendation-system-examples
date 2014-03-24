function [params,beta,time_params,nll] = max_arbp_obj_time(X,Y,adj,varargin)
%
% % example
% X = randn(200,40);
% Y = randn(100,20);
% adj = sprand(200,100,.1);
% [reg,beta,nll] = max_arbp_obj(X,Y,adj,'regr',lin_reg('ridge',1),'basis',@(x) exp(.1*x(:,1:2)));

% Tx=X(1,:);
% Ty=Y(1,:);

%% parameters
regr = lin_reg('ridge',.001);
maxiter = 100;
tol = 1e-4;
params = [];
basis = [];
ridge = eps;
learn_beta=1;
beta = [];
isstart = 1:size(Y,1);

SIMPLE_ARGS_DEF

%% time parameters
sel = find(~isstart);    
time_params = learn(exp_law,Y(sel,1)-Y(sel-1,1));  
%Y = Y(:,2:end);%remove the time output

%% initializations
[np,nc] = size(adj);
[I,J] = find(adj);%parent/child pairs
XX = X(I,:);
YY = Y(J,:);
XX(:,1)=YY(:,1)-XX(:,1);
YY=YY(:,2:end);
N = size(XX,1);
X(I,1)=XX(:,1);
basis_Xt = basis(X)';
%basis_Xt = basis(XX)';
if any(basis_Xt(:)<0)
    warning('there are some negative bases')    
end

    
c = sum(basis_Xt,2) + ridge;

R = size(basis_Xt,1);

if size(X,1)~=size(adj,1) ||  size(Y,1)~=size(adj,2)
    error('invalid sizes for parameters X,Y and adj')   
end

if isempty(params)
    params = update(regr,[],XX,YY,ones(N,1));
end

if isempty(beta)
    beta = ones(R,1);
end

ppa = adj2parents(adj);%cell array of parent indices
[ppa_vec,lbl] = cell2lbl(ppa);
lbl_cell = lbl2cell(lbl);
Ls = cellfun('length',ppa);

%% iterations
nll = inf*ones(1,maxiter);
w = eps;
tau = eps*ones(size(basis_Xt));
for it=2:maxiter
    
    %% tau
%    old_nll=nll;nll=neg_loglik(w,XX,YY,params,beta,tau,basis_Xt,regr,ridge,I);fprintf('%d) -ll = %4.2f(%3.1f)\n',1,nll,old_nll-nll);
%     if old_nll-nll<0
%         ddd
%     end
    lprob0 = logproba(regr,params,XX,YY);
%     temp=basis_Xt.*(beta*ones(1,N));
    temp=basis_Xt.*(beta*ones(1,np));
    tau = normalizedim(temp,1);
    
    %% w
%    old_nll=nll;nll=neg_loglik(w,XX,YY,params,beta,tau,basis_Xt,regr,ridge,I);fprintf('%d) -ll = %4.2f(%3.1f)\n',2,nll,old_nll-nll);
    lambda = max(1e-300,sum(temp,1)');
    [w,lp] = softmax_vec(lprob0+log(lambda(I)),lbl_cell,lbl);%faster with lbl pre-computed
     nll(it) = -sum(lp(Ls>0)) - logproba_prior(regr,params) + ridge*sum(beta) + sum(lambda);  
     fprintf('%d) -ll = %4.2f(%3.1f)\n',it-1,nll(it),nll(it-1)-nll(it));
     if abs(nll(it-1)-nll(it))<tol*abs(nll(it))
         break;
     end
    
    %% beta
%    old_nll=nll;nll=neg_loglik(w,XX,YY,params,beta,tau,basis_Xt,regr,ridge,I);fprintf('%d) -ll = %4.2f(%3.1f)\n',3,nll,old_nll-nll);
%     sw = accumarray(I,w,[N,1]); 
    sw = accumarray(I,w,[np,1]);
    if learn_beta
        beta =(tau*sw)./c;
    end
    
    %% params
%    old_nll=nll;nll=neg_loglik(w,XX,YY,params,beta,tau,basis_Xt,regr,ridge,I);fprintf('%d) -ll = %4.2f(%3.1f)\n',4,nll,old_nll-nll);
%    params = learn(regr,XX,YY,w);<=>params = update(regr,[],XX,YY,w);
    % TO DEBUG:
    % save bad_for_mlr.mat regr params XX YY w
    % load bad_for_mlr
    % nll=negloglik(regr,params,XX,YY,w);params2 = update(regr,params,XX,YY,w);nll2=negloglik(regr,params2,XX,YY,w);[nll-nll2]
    
    params = update(regr,params,XX,YY,w);
end




function nll = neg_loglik(w,XX,YY,params,beta,tau,basis_Xt,regr,ridge,I)
[R,n] = size(tau);
v1 = sum(w.*logproba(regr,params,XX,YY));
v2 = -sum(w.*log(w));
wpoint = accumarray(I,w,[n 1]);
v3 = sum(log(beta).*(tau*wpoint));
v4 = sum(sum(basis_Xt.*(beta*ones(1,n))));
v5 = wpoint'*sum(tau.*log(max(1e-300,basis_Xt)),1)';
v6 = -sum(sum(tau.*log(max(tau,1e-300)),1).*wpoint');
nll = -v1 - v2 - v3 + v4 - v5 - v6 - logproba_prior(regr,params) + ridge*sum(beta);
if ~isreal(nll)
    warning(' ')
end

