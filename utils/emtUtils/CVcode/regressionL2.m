function [model] = regressionL2(X,y,options)
% w = regressionL2(X,y)
%
% Returns Least Squares Solution:
%  min_w sum (Xw - y).^2
% by solving the Normal Equations
%
% options:
%   lambdaL2 = x (use L2 regularizer with non-negative strength x)

if nargin < 3
    options = [];
end

[verbose,lambdaL2,lambdaL1] = myProcessOptions(options,'verbose',1,'lambdaL2',0,'lambdaL1',0);
[nInstances,nVars] = size(X);

if lambdaL2 == 0 && lambdaL1 == 0
    % Least Squares
    w = (X'*X)\(X'*y);
elseif lambdaL1 == 0
    % Least Squares w/ L2 Prior (Ridge Regression)
    w = (X'*X + lambdaL2*eye(nVars))\(X'*y);
elseif lambdaL2 == 0
    % Least Squares w/ L1 Prior (Lasso)
    w = LassoShooting(X,y,lambdaL1,'verbose',verbose);
else
    % Least Squares w/ L2 and L1 Prior (Elastic Net)
    params.verbose = verbose*2;
    lambdaVect = repmat(lambdaL1,[nVars 1]);
    w = L1GeneralProjection(@penalizedL2,zeros(nVars,1),lambdaVect,params,@GaussianLoss,lambdaL2,X,y);
end

model.weights = w;
model.predictFunc = @(model,X)X*model.weights;
model.errL1Func = @(model,X,y)sum(abs(X*model.weights-y));
model.errL2Func = @(model,X,y)sum((X*model.weights-y).^2);
model.lossFunc = @(model,X,y)sum((X*model.weights-y).^2);
