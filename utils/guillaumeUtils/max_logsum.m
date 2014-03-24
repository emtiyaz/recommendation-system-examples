function [theta,objective] = max_logsum(X,nu,varargin)
% max_logsum - maximizes the sum of logarithms
%
% % example:
% %---------
% X=sparse(1:100,ceil(rand(1,100)*3),ones(1,100),100,3);
% X=normalizedim(sprand(100,3,.1),2);nu = rand(100,3);
% [theta,fopt] = max_logsum(X,nu);
% epsi=1;
% for i=1:1000
%   t2 = normalizedim(theta.*(1-epsi+rand(size(theta))*2*epsi),2);
%   f2 = logsum(X,nu,t2);dif(i) = f2-fopt;
%   fprintf('%4.3f\n',dif(i));
% end
% hist(dif,30);

[n,K] = size(X);
epsi=1;
objective = inf;
theta = normalizedim(ones(K,K),2); %intitialize theta
method = 'EM';
maxiter = 20;
verbose = 1;
prior_theta = .1;
tol = 1e-5;

SIMPLE_ARGS_DEF


switch(method)
    case 'random'
        for i=1:10000
            thetat =  normalizedim(theta.*(1-epsi+rand(size(theta))*2*epsi),2);
            ft = logsum(X,nu,thetat);
            if ft<objective
                theta = thetat;
                objective = ft;
            end            
        end
        
    case 'EM'
        objective = -inf;
        for it=1:maxiter

            objective_old = objective;
            %objective=logsum()+penalty term;
            objective = logsum_obj(X,nu,theta) + sum(sum(prior_theta.*log(theta)));
            delta=-objective_old+objective;
            if verbose
                %fprintf('%f\t',sum(sum(prior_theta.*log(theta))))
                fprintf('%d) nll=%4.3f(%3.2f)\n',it,objective,delta);
            end
            if abs(delta)<tol*(1-objective)
                %it debug
                break
            end
            
            %iterations
            
            W = normalizedim((X*theta').*nu,2);
            theta=normalizedim(W'*X+prior_theta,2);            

        end
    otherwise
        warning('unknown method')
end
