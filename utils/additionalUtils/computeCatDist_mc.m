function [postDist, Y1] = computeCatDist_mc(method, M, V, W, bias, K)
% Dec 5, 2011: bias included

  N = 10000;
  L = size(M,1);
  D = length(K);

  V = 0.5*(V + V');
  z = mvnrnd(M(:)', V, N);
  for d = 1:D
    idx = sum(K(1:d-1))+1:sum(K(1:d));
    eta_d = bsxfun(@plus, z*W(idx,:)',bias(idx)');
    switch method
    case {'log','boh','multiLogitLGGM','multiLogitFA','binaryLogitLGGM','binaryLogitFA'}
      lse = logsumexp(eta_d,2);
      prob = exp(bsxfun(@minus, eta_d, lse));
      prob = bsxfun(@times, prob, 1./sum(prob,2));
    case {'pw','stickLogitLGGM','stickLogitFA'}
      sigm = sigmoid(eta_d);
      sigm = sigm';
      prob = zeros(size(sigm));
      for m = 1:K(d)-1
        prob(m,:) = prod([(1-sigm(1:m-1,:)); sigm(m,:)], 1);
      end
      prob(end,:) = prod(1-sigm(1:K(d)-1,:), 1);
      prob = prob';
    otherwise
      error('no such method')
    end
    %Y(:,idx) = mnrnd(ones(N,1),prob);
    Y1 = mnrnd(ones(N,1),prob);
    Y(:,d) = sum(bsxfun(@times, Y1, [1:K(d)]),2); 
  end
  nidx = find(~isnan(sum(Y,2)));
  Y = Y(nidx,:);
  [Yunique,I,J]  = unique(Y,'rows'); 
  postDist = hist(J,1:size(Yunique,1))/size(Y,1);
  Y1 = sparse(encodeDataOneOfM(Yunique', K, 'M'));% one of K encoding

