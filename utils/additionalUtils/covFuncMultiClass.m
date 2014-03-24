function K = covFuncMultiClass(covfunc, nClass, hyp, x, z, class, i)
% covariance functions for multiclass GP. covfunc is an array of functions as
% listed in 'covFunctions.m'. nClass is number of classes. Rest of the elements
% are similar to the other GP covariance functions. K is a sparse matrix whose
% size depends on the mode.
% This functions currently works for the following covfuncs
%        {covConst,...
%        covLIN,...       
%        covLINone,...    
%        covNNone,...     
%        covNoise,...     
%        covPeriodic,...  
%        covRQiso,...     
%        covSEiso,...     
%        covSEisoU};    
% 
% %  s = covFuncMultiClass(covfunc, nClass) returns the # of hyperparameters required.
%  e.g. s = covFuncMultiClass({@covSEiso}, 3)
%
%  s = covFuncMultiClass(covfunc, nClass, hyp) is same as above.
%
% Other usage are similar to functions in covFunc.m
%  K = covFuncMultiClass(covfunc, nClass, hyp, x) 
%  K = covFuncMultiClass(covfunc, nClass, hyp, x, []) 
%  K = covFuncMultiClass(covfunc, nClass, hyp, x, 'diag') 
%  K = covFuncMultiClass(covfunc, nClass, hyp, x, xs) 
%  K = covFuncMultiClass(covfunc, nClass, hyp, x, [], i) 
%  K = covFuncMultiClass(covfunc, nClass, hyp, x, 'diag', i) 
%  K = covFuncMultiClass(covfunc, nClass, hyp, x, xs, i) 
%
% Written by Emtiyaz, CS, UBC
% Modified on Nov. 27, 2011

% some checks
if nargin >1
  M = nClass-1;
  if nargin >2
    if size(hyp,1) ~= M;
      error('hyp should be a matrix of size (nClass-1)x %s\n', feval(covfunc{:}));
    end
  end
end
if nargin <2
  error('not enough input');
  return
end

% determine mode
if nargin > 4
  if strcmp(z,'diag');
    mode = 1;
  elseif isempty(z)
    mode = 2;
  else 
    mode = 3;
  end 
end
if nargin < 4 
% return the number of hyperparameters
  K = M*str2num(feval(covfunc{:}));

elseif nargin == 4
% return K(x,x)
  N = size(x,1);
  K = zeros(N*M,N*M);
  for k = 1:M
    idx = [k:M:(N-1)*M+k];
    K(idx,idx) = feval(covfunc{:}, hyp(k,:), x);
  end
  K = sparse(K);

elseif nargin == 5
% return K(x,xs) or diag(K(x,x))
  N = size(x,1);
  switch mode
  case 1 % diag(K(x,x))
    K = zeros(N*M,1);
    for k = 1:M
      idx = [k:M:(N-1)*M+k];
      K(idx,1) = feval(covfunc{:}, hyp(k,:), x, z);
    end
  case 2 % K(x,x)
    K = zeros(N*M,N*M);
    for k = 1:M
      idx = [k:M:(N-1)*M+k];
      K(idx,idx) = feval(covfunc{:}, hyp(k,:), x, z);
    end
  case 3 % diag(K(x,xs))
    Ns = size(z,1);
    K = zeros(N*M,Ns*M);
    for k = 1:M
      idx1 = [k:M:(N-1)*M+k];
      idx2 = [k:M:(Ns-1)*M+k];
      K(idx1,idx2) = feval(covfunc{:}, hyp(k,:), x, z);
    end
  otherwise 
    error('no such mode')
  end
  K = sparse(K);

elseif nargin > 5
% return gradient
  N = size(x,1);
  idx = [class:M:(N-1)*M+class];
  switch mode
  case 1 % gradient of diag(K(x,x))
    K = zeros(N*M,1);
    K(idx,1) = feval(covfunc{:}, hyp(class,:), x, z, i);
  case 2 % gradient of K(x,x)
    K = zeros(N*M,N*M);
    K(idx,idx) = feval(covfunc{:}, hyp(class,:), x, z, i);
  case 3 % gradient of K(x,xs)
    Ns = size(z,1);
    K = zeros(N*M,Ns*M);
    idx1 = [class:M:(Ns-1)*M+class];
    K(idx,idx1) = feval(covfunc{:}, hyp(class,:), x, z, i);
  otherwise 
    error('no such mode')
  end
  K = sparse(K);

else
  error('too many output arguments');
end

