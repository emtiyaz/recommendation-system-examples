function logP = myMultNormPdf(X, mean_, covMat)
% x is dimxN where N is the number of samples and dim is dimensiion

  [dim, N] = size(X);

  mean_ = mean_(:);
  diff = chol(covMat)*(X - repmat(mean_,1,N));
  logP = -0.5*(sum(diff.^2,1));
  logP = logP - 0.5*logdet(2*pi*covMat);

