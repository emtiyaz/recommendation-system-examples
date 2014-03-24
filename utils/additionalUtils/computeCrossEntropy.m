function score = computeCrossEntropy(phat, y, nClass)

  % remove any zeros from phat
  phat = max(phat, eps);
  phat = bsxfun(@times, phat, 1./sum(phat,1));

  N = length(y);
  y = encodeDataOneOfM(y, nClass*ones(N,1), 'M');
  y = reshape(y, nClass, N);
  score = -mean(sum(log2(phat).*y,1));
