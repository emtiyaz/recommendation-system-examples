function pass = checkIfSumsToOne(A, dim)
% Checks if imput matrix sums to one or not. Only works for two dimensional
% matrices. 
% Written by Emtiyaz, XRCE, Grenoble, France 
% modified on Aug. 20, 2009

  tol = 1e-10;
  a = sum(A,dim);
  pass = 1;
  for n = 1:length(a)
    if a(n)>1+tol | a(n)<1-tol
      error('This matrix does not sum to one');
      pass = 0;
    end
  end


