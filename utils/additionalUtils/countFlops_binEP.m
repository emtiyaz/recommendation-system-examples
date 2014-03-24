function nflops = countFlops_binEP(N,D,nSweeps)

flops(0);

%% K = feval(cov{:}, hyp.cov, x); 
%addflops(countFlops_covSEiso(N,D));
%nlZ0 = -sum(feval(lik, hyp.lik, y, m, diag(K), inf));
%addflops(countFlops_likLogistic(N));
%while (abs(nlz-nlz_old) > tol && sweep < max_sweep) || sweep<min_sweep
addflops(nSweeps*countFlops_epSiteUpdate(N));
%sW = sqrt(ttau); 
% alpha = tnu-sW.*solve_chol(L,sW.*(K*tnu));
%addflops(N*flops_sqrt);
%addflops(flops_mul(N,N,1) + N + flops_chol(N) + 2*flops_solve_tri(N,N,1));

nflops = flops;
return;

function nflops = countFlops_epSiteUpdate(N);

flops(0);
% for i = 1:n%randperm(n)
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
% [Sigma, mu, nlZ, L] = epComputeParams(K, y, ttau, tnu, lik, hyp, m, inf);
nflops = N*flops;
nflops = nflops + countFlops_epComputeParams(N);

return;

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

%lZ = feval(lik, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf);
%addflops(countFlops_likLogistic(N));
%nlZ = sum(log(diag(L))) -sum(lZ) -tnu'*Sigma*tnu/2  ...
  %  -(nu_n-m.*tau_n)'*((ttau./tau_n.*(nu_n-m.*tau_n)-2*tnu)./(ttau+tau_n))/2 ...
  %  +sum(tnu.^2./(tau_n+ttau))/2-sum(log(1+ttau./tau_n))/2;
%addflops(2*N - N + flops_mul(N,N,1) + flops_mul(1,N,1));
%addflops(2*N + 3*N + flops_mul(1,N,1));
%addflops(3*N + N);
 
nflops = flops;
return;

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
      return;
    
function nflops = countFlops_covSEiso(N,D)

flops(0);
%ell = exp(hyp(1)); sf2 = exp(2*hyp(2));
addflops(2*flops_exp)
%K = sq_dist(x'/ell,z'/ell);
addflops(N*D + flops_mul(N,D,N));
%K = sf2*exp(-K/2);
addflops(N^2 + N^2*flops_exp);

nflops = flops;
return;

