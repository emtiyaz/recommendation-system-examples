function [theta,objective] = max_arbp_obj_disc(X,nu,adj,varargin)
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
[nhist] = size(nu,1);
objective_old = -inf;
theta = []; %intitialize theta
maxiter = 20;
verbose = 1;
tol = 1e-5;
classifier = multinologreg_cl('ridge',1e-3);

SIMPLE_ARGS_DEF

if isempty(theta)
    theta = learn(classifier,zeros(size(nu)),normalizedim(ones(n,K),2));
    %theta = normalizedim(ones(K,K),2);
end
%save('test_max_arbp_obj','X','nu','adj','varargin');

Xp=X';
pa = adj2parents(adj);%cell array of parent indices
%child = adj2parents(adj');%cell array of child indices
Ls = cellfun('length',pa);
L = max(Ls);%maximal nb of parents

[nul,Xcat] = max(Xp,[],1);%categories

% I = zeros(L,n);
% for i=1:n
%     cL = length(pa{i})
%     I(1:cL,i) = pa{i};
%     J(1:cL,i) = 1:cL;
% end
% sbs = sub2ind([n,K],I',J');

%% EM algorithm
for it=1:maxiter
    %E step
    lp = logproba(classifier,theta,nu)';%nhist*K matrix
%    W0 = lp(sbs);%n*L matrix
    W0 = -inf*ones(n,L);
    for i=1:n
        if ~isempty(pa{i})
            W0(i,1:length(pa{i})) = Xp(:,i)'*lp(:,pa{i});
        end
    end 
    [W,ll] = softmax(W0',1);
    
        %objective=logsum()+penalty term;
        objective = sum(ll(Ls>0)) + logproba_prior(classifier,theta);
        delta = objective-objective_old;
        if verbose
            fprintf('%d) nll=%4.3f(%3.2f)\n',it,objective,delta);            
        end
        if abs(delta)<tol*(1+abs(objective))
            break
        end
        objective_old = objective;
    % M step
%    Y = accumarray(asgn,W,n,K);
    Y  = zeros(nhist,K);
    for i=1:n
        if ~isempty(pa{i})
           Y(pa{i},Xcat(i)) = Y(pa{i},Xcat(i)) + W(1:length(pa{i}),i);
        end
    end    
    theta = learn(classifier,nu,Y);
end



