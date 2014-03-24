function out = scalar2Vec(in,L)
% converts scalar to a vector of length L

  if isscalar(in);
    out = ones(L,1)*in;
  else 
    out = in;
  end

