function mu = meanFunctionsMultiClass(meanfunc, nClass, hyp, x, class, i)
% Mean functions for multiclass GP. meanfunc is an array of functions as
% listed in 'meanFunctions.m'. nClass is number of classes. Rest of the elements
% are similar to the other GP covariance functions. 
% This functions currently works for the following meanfuncs
%       {@meanZero,...   
%         @meanOne,...    
%         @meanConst,...  
%         @meanLinear}; 
% 
% %  s = meanFuncMultiClass(meanfunc, nClass) returns the # of hyperparameters required.
%  e.g. s = meanFuncMultiClass({@meanZero}, 3)
%
%  s = meanFuncMultiClass(meanfunc, nClass, hyp) is same as above.
%
% Other usage are similar to functions in meanFunc.m
%  mu = meanFuncMultiClass(covfunc, nClass, hyp, x) 
%  mu = meanFuncMultiClass(covfunc, nClass, hyp, x, i) 
%
% Written by Emtiyaz, CS, UBC
% Modified on Nov. 27, 2011

% some checks
if nargin >1
  M = nClass-1;
end
if nargin <2
  error('not enough input');
  return
end

if nargin < 4 
% return the number of hyperparameters
  mu = M*str2num(feval(meanfunc{:}));

elseif nargin == 4
% return mu(x)
  if isempty(hyp); hyp = zeros(M,0); end;
  N = size(x,1);
  mu = zeros(N*M,1);
  for k = 1:M
    idx = [k:M:(N-1)*M+k];
    mu(idx) = feval(meanfunc{:}, hyp(k,:)', x);
  end

elseif nargin > 4
% return gradient
  if isempty(hyp); hyp = zeros(M,0); end;
  N = size(x,1);
  idx = [class:M:(N-1)*M+class];
  mu = zeros(N*M,1);
  mu(idx) = feval(meanfunc{:}, hyp(class,:)', x, i);

else
  error('too many output arguments');
end

