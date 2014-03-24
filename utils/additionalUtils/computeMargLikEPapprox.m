function [nlZ] = computeMargLikEPapprox(post, hyp, mean, cov, lik, x, y)
% computes EP approx to margLik 
% post is posterior distribution, rest is same input as other gpml functions
% Also, y is -1 or 1, 
% For usage, see testInferBernLogitGP_pw.m
% Written by Emtiyaz, CS, UBC, 
% Modified on Sep 27, 2011

  lambda = (post.sW).^2;
  alpha = post.alpha;
  Sigma = post.covMat;
  L = post.L;
  mean_ = feval(mean{:}, hyp.mean, x);
  covMat = feval(cov{:}, hyp.cov, x); 
  m = mean_;

  % compute mu_i, then sigma_\i and mu_\i 
  % see Kuss and Rasmussen, 2004 Eq. 24
  mu = (diag(1./lambda) + covMat)*alpha;
  
  % Fix this
  v_i = 1./(1./(diag(post.covMat)) - lambda);
  
  % S:debug
  a = find(v_i < 0);
  if (length(a) > 0)
      a
      %error('Too many neg v_i');
  end;
  v_i(v_i <0) = eps;
  
  
  m_i = post.mean./diag(post.covMat) - mu.*lambda;
  m_i = v_i.*m_i;

  % This part is copied from infEP function
  ttau = lambda;
  tnu = lambda.*mu;
  nu_n = m_i./v_i;
  tau_n = 1./v_i;
  % compute log(Z_i)
  lZ = feval(lik, hyp.lik, y, nu_n./tau_n, 1./tau_n, 'infEP');
  % compute nlZ
  nlZ = sum(log(diag(L))) -sum(lZ) -tnu'*Sigma*tnu/2  ...
    -(nu_n-m.*tau_n)'*((ttau./tau_n.*(nu_n-m.*tau_n)-2*tnu)./(ttau+tau_n))/2 ...
    +sum(tnu.^2./(tau_n+ttau))/2-sum(log(1+ttau./tau_n))/2;

