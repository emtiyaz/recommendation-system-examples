function lvp = logsumexpbnd_optim_improved(mu,var,lvp,maxiter,method)
% logsumexpbnd_optim_improved - optimize the improved log-sum-exp bound 
%


[n,K] = size(mu);

if  nargin<3|| isempty(lvp)
    lvp = [];
    T = eye(K);
else
    T = eye(K)-ones(K,1)*lvp.beta';
end

if nargin<4 || isempty(maxiter);
    maxiter = 100;
end

if nargin<5
    method = 0;
end



% if ~isempty(lvp)
%     lvp = rmfield(lvp,'beta');
% end


lvp_old = lvp;

lvp = logsumexpbnd_optim((T*mu')',diag(T*var*T')',lvp,maxiter,method);

%check convergence
% lvp_old.method = 3;
% lvp.method = 3;
% check_update(lvp_old,lvp,mu,var)
% lvp_old = lvp;

% lvpNoBeta = rmfield(lvp,'beta');
lvpNoBeta = lvp;lvpNoBeta.beta=zeros(K,1);
[A,b,c] = logsumexpbnd_poly(lvpNoBeta,1); %#ok<NASGU>
E2 = var + mu'*mu;
lvp.beta = -1./2./(sum(sum(A))) * pinv(E2)*(mu' - mu'*sum(b) - 2*E2*A*ones(K,1));

% %check convergence
% lvp_old.method = 3;
% lvp.method = 3;
% check_update(lvp_old,lvp,mu,var)



















%% check validity of updates
function check_update(lvp_old,lvp,mu,var)

if ~isempty(lvp_old) && isfield(lvp_old,'chi')
    [a_old,b_old,c_old] = logsumexpbnd_poly(lvp_old,1);        
    E_old = mu*a_old*mu' + sum(sum(var.*a_old)) + mu*b_old' + c_old;

    % % equivalent
%     beta = lvp_old.beta;
%     K = length(beta);
%     T = eye(K)-ones(K,1)*beta';
%     lvp_old.beta = lvp_old.beta*0;
%     [a_old,b_old,c_old] = logsumexpbnd_poly(lvp_old,1);        
%     E_old2 = (mu*T')*a_old*(T*mu') + sum(sum((T*var*T').*a_old)) + mu*T'*b_old' + c_old + beta'*mu'
else
    E_old = inf;
end
if ~isempty(lvp) && isfield(lvp,'chi')
    [a,b,c] = logsumexpbnd_poly(lvp,1);    
    E = mu*a*mu' + sum(sum(var.*a)) + mu*b' + c;   
    fprintf('%20s: %f (%5f)\n','value',E,E_old-E);
    if E_old-E<-1e-2 
        save temp.mat lvp_old lvp mu var
        lvp_old.beta
        lvp.beta
        error('bound increase!')
    end
end

