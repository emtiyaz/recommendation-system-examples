function nflops = countFlops(path, algo, params)

switch algo
    case 'girolami_vb'
        D = params.D;
        nClass = params.nClass;
        nIters = path;
        
        iter_flops = countFlops_girolamviVBMP(D, nClass, nIters);
        % returns cumulative sum
        nflops = cumsum(iter_flops);
        nflops = [0; nflops(:)];
        
    case 'ep'
        D = params.D;
        nflops = countFlops_epSiteUpdate(D);
        
    case 'naive'
        D = params.D;
        R = params.R;
        fevals = path.funcCount;
        nIters = length(fevals);
        nflops = zeros(nIters,1);
        func_flops = countFlops_funObj_varGaussApproxLGGM_noReparam(D,R);
        % returns cumulative sum
        nflops = fevals*func_flops;
        
    case 'coordAscentReparam'
        D = params.D;
        R = params.R;
        fevals = path.fevals;
        fevals_alpha = fevals.alpha;
        fevals_lambda = fevals.lambda;
        func_alpha_flops = countFlops_funObjM_reparam(D,R);
        func_lambda_flops = countFlops_funObjV_reparam(D,R);
        nflops = func_alpha_flops*fevals_alpha + func_lambda_flops*fevals_lambda;
        nflops = cumsum(nflops);
        
    case {'reparam_orig','reparam'}
        D = params.D;
        R = params.R;
        fevals = path.funcCount;
        nIters = length(fevals);
        nflops = zeros(nIters,1);
        func_flops = countFlops_funObj_varGaussApproxLGGM_reparam_orig(D,R);
        % returns cumulative sum
        nflops = fevals*func_flops;
        
    case 'fixedPointMultiClass'
        D = params.D;
        R = params.R;
        nClass = params.nClass;
        ind = params.ind;
        
        fevals = path.fevals;
        fevalsFixPt = fevals(:,2);
        fevalsM = fevals(:,1);
        nIters = size(fevals,1);
        
        nflops = zeros(nIters,1);
        % update V
        nflops_updateV = countFlops_fixPt_updateV_multiClass(D,nClass);
        nflops = nflops + nflops_updateV;
        % fixed point solver
        nflops_fixPt = countFlops_emFixPt_solveFixedPoint(R);
        nflops = nflops + fevalsFixPt*nflops_fixPt;
        % update mean
        func_flopsM = countFlops_fixPt_funcMean_multiClass(D,nClass,R,ind);
        nflops = nflops + fevalsM*func_flopsM;
        % returns cumulative sum
        nflops = cumsum(nflops);
        nflops = [0; nflops(:)];
        
    case {'fixedPoint','fixedPointv1'}
        D = params.D;
        R = params.R;
        fevals = path.fevals;
        fevalsFixPt = fevals(:,2);
        fevalsM = fevals(:,1);
        nIters = size(fevals,1);
        
        nflops = zeros(nIters,1);
        % update V
        if strcmp(algo,'fixedPointv1')
            nflops_updateV = countFlops_emFixPt_updateV_v1(D);
        else
            nflops_updateV = countFlops_emFixPt_updateV(D);
        end;
        
        nflops = nflops + nflops_updateV;
        % fixed point solver
        nflops_fixPt = countFlops_emFixPt_solveFixedPoint(R);
        nflops = nflops + fevalsFixPt*nflops_fixPt;
        % update mean
        func_flopsM = countFlops_emFixPt_funcMean(D,R);
        nflops = nflops + fevalsM*func_flopsM;
        % returns cumulative sum
        nflops = cumsum(nflops);
        nflops = [0; nflops(:)];
        
    case 'emFixPt'
        [D,N] = size(Y);
        Dobs = sum(Y ~= 0, 1); Dobs = Dobs(:);
        Nobs = sum(Y ~= 0, 2); Nobs = Nobs(:);
        
        % 1st row contains nFevals for funcMean
        % 2nd row contains nFevals for fixed point solver, and so on...
        feval_Estep = path.Estep;
        feval_Mstep = path.Mstep;
        nMVupdate = 4;
        
        nIters = size(feval_Mstep,1);
        R = 20;% nPieces
        nflops_fixPt = countFlops_emFixPt_solveFixedPoint(R);
        nflops = zeros(nIters,1);
        for i = 1:nIters
            for n = 1:N
                Dn = Dobs(n);
                % Omega = inv(covMat(obs,obs))
                nflops(i) = nflops(i) + flops_inv(Dn);
                % update V
                nflops_updateV = countFlops_emFixPt_updateV(Dn);
                nflops(i) = nflops(i) + nMVupdate*nflops_updateV;
                % cost of fixed point solver
                nflops(i) = nflops(i) + feval_Estep(2*(i-1)+2,n)*nflops_fixPt;
                % compute lambda
                nflops(i) = nflops(i) + 173*R + Dn;
                % update m
                func_flopsM = countFlops_emFixPt_funcMean(Dn,R);
                nflops(i) = nflops(i) + feval_Estep(2*(i-1)+1,n)*func_flopsM;
                % comput ss
                nflops(i) = nflops(i) + countFlops_emReparam_ss(Dn);
            end
            % M-step
            nflops(i) = nflops(i) + countFlops_emReparam_mstep(D);
        end
        
    case {'emReparam','emReparamApprox'}
        feval_Estep = path.Estep;
        nIters = size(feval_Estep,1);
        nflops = zeros(nIters,1);
        R = 20;% nPieces
        for i = 1:nIters
            for n = 1:N
                Dn = Dobs(n);
                % funObj counts
                func_flops = countFlops_funObj_varGaussApproxLGGM_reparam_orig(Dn,R);
                nflops(i) = nflops(i) + feval_Estep(i,n)*func_flops;
                % comput ss
                nflops(i) = nflops(i) + countFlops_emReparam_ss(Dn);
            end
            % M-step
            nflops(i) = nflops(i) + countFlops_emReparam_mstep(D);
        end
    otherwise
        error('no such method');
end

function nflops = countFlops_funObjV_reparam(D,R);
flops(0);
%lambda = exp(x);
addflops(flops_exp*D)
%slambda = sqrt(lambda);
addflops(D*flops_sqrt)
%B = eye(N) + (slambda*slambda').*covMat;
addflops(D + flops_mul(D,1,D) + (D^2))
%U = chol(B);
addflops(flops_chol(D))
%Uinv = U\eye(N);
addflops(flops_solve_tri(D,D,D))
%A = diag(1./slambda)*Uinv*Uinv'*diag(slambda);
addflops(D^2 + flops_mul(D,D,D))
%V = A*covMat;
addflops(flops_mul(D,D,D))
%kl = -0.5*(2*sum(log(diag(U))) + trace(A));
addflops(1 + D + D*flops_log + D)
%[fb, gm, gv] = funObj_pw(M, diag(V), bound);
addflops(D*R*172);
%f = kl - sum(fb);
addflops(1+D);
%g = (V.^2)*(-0.5*lambda + gv);
addflops(D^2 + D + 1 + flops_mul(D,D,1));
%g = lambda.*g;
addflops(D);
nflops = flops;

function nflops = countFlops_funObjM_reparam(D,R);
flops(0)
%M = mean_ + covMat*alpha;
addflops(D+flops_mul(D,D,1));
%kl = -0.5*(alpha'*covMat*alpha);
addflops(1+flops_mul(1,D,1));
%[fb, gm, gv] = funObj_pw(M, diag(V), bound);
addflops(D*R*172);
%f = kl + y'*M - sum(fb);
addflops(2 + D + D);
%g = covMat*(-alpha + y - gm);
addflops(2*D + flops_mul(D,D,1));
nflops = flops;

function nflops = countFlops_emFixPt_solveFixedPoint(R)
flops(0)
%[fb, gmb, gvb] = funObj_pw_vec(m, x, bound, [1 0 1]);
addflops(R*172);
%f(i) = log(x) + (t-Omega)*x - 2*fb;
addflops(flops_log + 4);
%k = Omega + 2*gvb;
%x = 1/(k-t);
addflops(5);
nflops = flops;

function nflops = countFlops_fixPt_funcMean_multiClass(D,nClass,R,ind);
% flops for one funcMean evals

M = nClass-1;
DM = D*M;
flops(0);
%[fb, gmb, gvb] = funObj_pw_vec(m, v, bound, [1 1 0]);
addflops(sum(ind)*R*172);
%diff = (m-mu);
addflops(DM);
%gm = Omega*diff; Omega is block structured
addflops(M*flops_mul(D,D,1));
%f = 0.5*diff'*gm - y'*m + sum(ind.*fb);
addflops(2*flops_mul(1,DM,1) + sum(ind) + 3);
%g = gm - y + ind.*gmb;
addflops(2*DM);

nflops = flops;

function nflops = countFlops_emFixPt_funcMean(D,R);
% flops for one funcMean evals
flops(0);
%[fb, gmb, gvb] = funObj_pw_vec(m, v, bound, [1 1 0]);
addflops(D*R*172);
%diff = (m-mu);
addflops(D);
%gm = Omega*diff;
addflops(flops_mul(D,D,1));
%f = 0.5*diff'*gm - y'*m + sum(fb);
addflops(2*flops_mul(1,D,1) + D + 3);
%g = gm - y + gmb;
addflops(2*D);
nflops = flops;

function nflops = countFlops_fixPt_updateV_multiClass(D,nClass);

DM = D*(nClass-1);
% flops for one update of V (except the fixed point solution)
flops(0)
%t11 = v12/sqrt(v22);
addflops(DM-1)
%K11inv = V11 - t11*t11';
addflops(2*(DM-1)^2)
%t1 = K11inv*o12; since o12 is sparse
%t = o12'*t1;
addflops(flops_mul(D,D,1) + flops_mul(1,D,1))
%t1 = t1*sqrt(v22);
addflops(DM-1 + flops_sqrt)
%v12 = -t1*sqrt(v22);
addflops(DM-1)
%V(idx,idx) = K11inv + t1*t1';
addflops(2*(DM-1)^2)
% repeat for each dimension d and category k
nflops = DM*flops;

function nflops = countFlops_emFixPt_updateV(D);
% flops for one update of V (except the fixed point solution)
flops(0)
%t11 = v12/sqrt(v22);
addflops(D-1)
%K11inv = V11 - t11*t11';
addflops(2*(D-1)^2)
%t1 = K11inv*o12;
addflops(flops_mul(D-1,D-1,1))
%t = o12'*t1;
addflops(flops_mul(1,D-1,1))
%t1 = t1*sqrt(v22);
addflops(D-1 + flops_sqrt)
%v12 = -t1*sqrt(v22);
addflops(D-1)
%V(idx,idx) = K11inv + t1*t1';
addflops(2*(D-1)^2)
% repeat for each dimension d
nflops = D*flops;

function nflops = countFlops_emFixPt_updateV_v1(D);
% flops for one update of V (except the fixed point solution)
flops(0)
%t = Kdiag(i) - 1/v22;
addflops(2)
%V(idx,idx) = V(idx,idx) + ((v22-v22_old)/(v22_old)^2)*(v12*v12');
addflops(2*(D-1)^2 + (D-1))
%V(idx,i) = -(v22/v22_old)*v12;
addflops(D-1)
% repeat for each dimension d
nflops = D*flops;

function nflops = countFlops_emReparam_mstep(D)

flops(0)
%delta = (covMat*ssM)/N;
addflops(flops_mul(D,D,1) + D);
%mean_ = mean_ + delta;
addflops(D);
%sigma_emp = covMat + covMat*ssC*covMat/N - delta*delta';
addflops(2*D^2 + flops_mul(D,D,D) + flops_mul(D,1,D));
nflops = flops;

function nflops = countFlops_emReparam_ss(Dn)
%
flops(0)
%lambda = exp(v(Dn+1:end));
addflops(flops_exp*Dn)
%ssM(obs) = ssM(obs) + alpha;
addflops(Dn)
%slambda = sqrt(lambda);
addflops(Dn*flops_sqrt)
%U = chol(eye(Dn) + (slambda*slambda').*covMat(obs,obs));
addflops(flops_mul(Dn,1,Dn) + (Dn^2)*flops_mul(1,1,1))
addflops(flops_chol(Dn))
%Uinv = U\eye(Dn);
addflops(flops_solve_tri(Dn,Dn,Dn))
%T = bsxfun(@times, Uinv', slambda(:)');
addflops(Dn^2*flops_mul(1,1,1))
%invB = T'*T;
addflops(flops_mul(Dn,Dn,Dn))
%ssC(obs,obs) = ssC(obs,obs) + alpha*alpha' - invB;
addflops(flops_mul(Dn,1,Dn))
addflops(3*Dn^2)
nflops = flops;
return;

function nflops = countFlops_funObj_varGaussApproxLGGM_reparam_orig(Dn,R)

%%%%%%%%%%%
% funObj
%%%%%%%%%%%
% funObj_varGaussReparam_orig
flops(0);
%lambda = exp(x(N+1:end));
addflops(flops_exp*Dn)
%m = mean_ + covMat*alpha;
addflops(flops_mul(Dn,Dn,1) + Dn);
%slambda = sqrt(lambda);
addflops(flops_sqrt*Dn);
%B = eye(N) + (slambda*slambda').*covMat;
addflops(flops_mul(Dn,1,Dn) + (Dn^2) + Dn)
%U = chol(B);
addflops(flops_chol(Dn))
%Uinv = U\eye(N);
addflops(flops_solve_tri(Dn,Dn,Dn))
%A = diag(1./slambda)*Uinv*Uinv'*diag(slambda);
addflops(flops_mul(Dn,Dn,Dn))
addflops(Dn^2)
%V = A*covMat;
addflops(flops_mul(Dn,Dn,Dn))
%kl = -0.5*(2*sum(log(diag(U))) + trace(A) + alpha'*covMat*alpha -N);
addflops(flops_log*Dn + 2*Dn + flops_mul(Dn,Dn,1) + flops_mul(1,Dn,1) + 4 + 2);
%[f, gm, gv] = funObj_pw(m, diag(V), bound);
addflops(Dn*R*219);
%f = kl + y'*m - sum(f);
addflops(flops_mul(1,Dn,1) + Dn + 2)
%galpha = covMat*(-alpha + y - gm);
addflops(flops_mul(Dn,Dn,1) + 2*Dn)
%glambda = (V.^2)*(-0.5*lambda + gv);
addflops(Dn^2 + flops_mul(Dn,Dn,1) + 2*Dn)
%glambda = lambda.*glambda;
addflops(Dn)
nflops = flops;


function nflops = countFlops_funObj_varGaussApproxLGGM_noReparam(D,R);

flops(0)
%[U,posdef] = cholproj(V);
addflops(flops_chol(D))
%logdetV = 2*sum(log(diag(U)));
addflops(flops_log*D +D+1)
%invV = solve_chol(U, eye(D));%U\(U'\eye(L));
addflops(flops_solve_tri(D,D,D))
%MMu = M-Mu;
addflops(D)
%invSigmaMMu = invSigma*MMu;
addflops(flops_mul(D,D,1))
%f = f - 0.5*(logdetSigma-logdetV + V(:)'*invSigma(:) + MMu'*invSigmaMMu - D);
addflops(1 + 1 + 4 + flops_mul(1,D^2,1) + flops_mul(1,D,1))
%gM = gM - invSigmaMMu;
addflops(D)
%gV = gV + 0.5*(invV-invSigma);
addflops(2*D^2)
%f = f + y'*mbar;
addflops(1+D)
%gM = gM + y;
addflops(D)
%[fb, gmb, gvb] = funObj_pw(mbar, vbar, bound);
addflops(D*R*219);
%f = f - sum(fb);
addflops(D+1);
%gM = gM - gmb;
addflops(D);
%gV = gV - diag(gvb); %gV = gV - W'*diag(gvb)*W;
addflops(D);
%gV = gV + triu(gV,1);
addflops(ceil(D*(D-1)/2));
nflops = flops;

function nflops = countFlops_epSiteUpdate(N);
% for EP site updates
% for i = 1:n%randperm(n)
flops(0);
% tau_ni = 1/Sigma(i,i)-ttau(i);      %  first find the cavity distribution ..
addflops(2);
% nu_ni = mu(i)/Sigma(i,i)+m(i)*tau_ni-tnu(i);    % .. params tau_ni and nu_ni
addflops(4);
%[lZ, dlZ, d2lZ] = feval(lik, hyp.lik, y(i), nu_ni/tau_ni, 1/tau_ni, inf);
addflops(countFlops_likLogistic(1));
% ttau(i) = -d2lZ  /(1+d2lZ/tau_ni);
addflops(3)
% ttau(i) = max(ttau(i),0); % enforce positivity i.e. lower bound ttau by zero
addflops(1)
% tnu(i)  = ( dlZ + (m(i)-nu_ni/tau_ni)*d2lZ )/(1+d2lZ/tau_ni);
addflops(7)
%ds2 = ttau(i) - ttau_old;                   % finally rank-1 update Sigma ..
addflops(1)
%Sigma = Sigma - ds2/(1+ds2*si(i))*si*si';          % takes 70% of total time
addflops((flops_mul(N,1,N) + 4 + N^2));
%mu = Sigma*tnu;                                        % .. and recompute mu
addflops(flops_mul(N,N,1));
%end
nflops = N*flops;
% [Sigma, mu, nlZ, L] = epComputeParams(K, y, ttau, tnu, lik, hyp, m, inf);
nflops = nflops + countFlops_epComputeParams(N);

function nflops = countFlops_epComputeParams(N)
flops(0);
%ssi = sqrt(ttau);
addflops(N*flops_sqrt)
%L = chol(eye(n)+ssi*ssi'.*K);
addflops(N + flops_mul(N,1,N) + N^2 + flops_chol(N));
%V = L'\(repmat(ssi,1,n).*K);
addflops(N^2 + flops_solve(N,N,N));
%Sigma = K - V'*V;
addflops(N^2 + flops_mul(N,N,N));
% mu = Sigma*tnu;
addflops(flops_mul(N,N,1));
%tau_n = 1./diag(Sigma)-ttau;
addflops(2*N);
%nu_n  = mu./diag(Sigma)-tnu+m.*tau_n;
addflops(4*N);

nflops = flops;

function nflops = countFlops_likLogistic(N)

flops(0);

% y = y.*ones(size(mu));
addflops(N);
% likLogistic(t) \approx 1/2 + \sum_{i=1}^5 (c_i/2) erf(lam_i/sqrt(2)t)
%lam = sqrt(2)*[0.44 0.41 0.40 0.39 0.36];    % approx coeffs lam_i and c_i
%c = [1.146480988574439e+02; -1.508871030070582e+03; 2.676085036831241e+03;
%    -1.356294962039222e+03;  7.543285642111850e+01                      ];
addflops(5);
% [lZc,dlZc,d2lZc] = likErf([], y*ones(1,5), mu*lam, s2*(lam.^2), inf);
addflops(N*5);
%lZ = log_expA_x(lZc,c);
%maxA = max(A,[],2); y = log(exp(A-maxA*ones(1,N))*x) + maxA;
addflops(N*5);
addflops(5*N + flops_mul(N,1,5) + 5*N*flops_exp + flops_mul(N,5,1) + N + N);
%dlZ  = expABz_expAx(lZc, c, dlZc, c.*lam');
% A = A-maxA*ones(1,N); y = ( (exp(A).*B)*z ) ./ ( exp(A)*x );
addflops(N*5);
addflops(N*5*flops_exp + N*5 + flops_mul(N,5,1) + flops_mul(N,5,1) + N*flops_div);
%d2lZ = expABz_expAx(lZc, c, dlZc.^2+d2lZc, c.*(lam.^2)') - dlZ.^2;
addflops(N*5);
addflops(N*5*flops_exp + N*5 + flops_mul(N,5,1) + flops_mul(N,5,1) + N*flops_div);
%val = abs(mu)-196/200*s2-4;
addflops(4);
%lam = 1./(1+exp(-10*val));
addflops(1 + N*flops_exp + N*flops_div);
%lZtail = min(s2/2-abs(mu),-.1);
addflops(N);
%dlZtail = -sign(mu);
addflops(N);
%id = y.*mu>0;
% lZtail(id) = log(1-exp(lZtail(id)));
addflops(2*N);
addflops(N*flops_exp + N);
%lZ  = (1-lam).* lZ + lam.* lZtail;
addflops(2*N);
%dlZ = (1-lam).*dlZ + lam.*dlZtail;
addflops(2*N);

nflops = flops;

