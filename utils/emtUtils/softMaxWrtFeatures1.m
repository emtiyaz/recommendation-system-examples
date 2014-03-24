function [nll,g,H] = softMaxWrtFeatures1(z,B,c)
% evaluates nll, grad, and hessian wrt to the feature vector.z of the following
% function:-log(y=k/B,z) = -beta_k^T*z + logsumexp(beta_j^T*z)
% B is the parameter vector of size nClass-1 x numOfFeatures 
% c is the observed class, scalar value between 1:nClass
% z is the feature vetor
% Written by Emtiyaz, CS, UBC
% modified on March 12, 2010

[M,D] = size(B);

B = [B; zeros(1,D)];
Bz = B*z;
sumexpBz = sum(exp(Bz));
nll = -Bz(c) + log(sumexpBz); 
if nargout > 1
  pk = exp(Bz)./sumexpBz;
  p = sum(bsxfun(@times, B', pk(:)'),2);
  g = -B(c,:)' + p;
  %g = -sum(B(y,:))' + N*tg./tSum;
  %BExpBz = B'.*repmat(expBz',D,1);
  %g = -sum(B(y,:))' + N*tg./tSum;
end

if nargout > 2
  H = B'*diag(pk)*B - p*p';
end

