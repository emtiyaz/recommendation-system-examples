function [theta,fmin] = max_arbp_obj_old(X,nu,adj,varargin)
% max_arbp_obj - maximizes the sum of logarithms
%
% % example:
% %---------
% 
% n=100;
% X = sparse(1:n,ceil(rand(1,n)*3),ones(1,n),n,3);%data
% nu = randn(n,3);%histories
% L=5;%order
% adj=full(sparse((1:n-L)'*ones(1,L)+ones(n-L,1)*(0:L-1),(1:n-L)'*ones(1,L)+L,ones(n-L,L),n,n));
% [theta,fopt] = max_arbp_obj(X,nu,adj);
%
%
% epsi=1;
% for i=1:1000
%   t2 = normalizedim(theta.*(1-epsi+rand(size(theta))*2*epsi),2);
%   f2 = logsum(X,nu,t2);dif(i) = f2-fopt;
%   fprintf('%4.3f\n',dif(i));
% end
% hist(dif,30);


%% default parameters
[n,K] = size(X);
[nhist,dim_inputs] = size(nu);
fold = inf;
theta = []; %intitialize theta
maxiter = 20;
verbose = 1;
tol = 1e-5;
regr = lin_reg('ridge',1);

SIMPLE_ARGS_DEF

if isempty(theta)
    theta = learn(regr,nu,normalizedim(rand(n,K),2));%initialize with random outputs
end
%save('test_max_arbp_obj','X','nu','adj','varargin');%if it does not work

Xp=X';

%% adjacency matrix precomputation
[I,J] = find(adj);%parent/child pairs
inds = sub2ind(size(adj),I,J);
pa = adj2parents(adj);%cell array of parent indices
%child = adj2parents(adj');%cell array of child indices
Ls = cellfun('length',pa);
L = max(Ls);%maximal nb of parents


%% EM algorithm
for it=1:maxiter

    %E step
    lp = logproba(regr,theta,nu(I,:),X(J,:))';
    
    lp_adj = sparse(I,J,lp,nhist,n);
    [W,ll] = softmax(lp_adj,1);

%    [Wvec,ll] = softmax_vec(lp(:),pa);
 
    %fmin=logsum()+penalty term;
    
    fmin = -sum(ll(Ls>0)) - logproba_prior(regr,theta);
    if isnan(fmin)
        warning('a')
    end
    
    delta=fold-fmin;
    if verbose
        fprintf('%d) nll=%4.3f(%3.2f)\n',it,fmin,delta);
    end
    if abs(delta)<tol*(1+abs(fmin))
        break
    end
    fold = fmin;

    % M step
%    W(inds)'
    theta = learn(regr,nu(I,:),X(J,:),W(inds));
%     theta = learn(regr,nu(I,:),X(J,:),Wvec);
end



