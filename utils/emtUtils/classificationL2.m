function model = classificationL2(X, y, options)
% w = classificationL2(X,y)
%
if nargin < 3
    options = [];
end

[verbose,lambdaL2,lambdaL1] = myProcessOptions(options,'verbose',1,'lambdaL2',0,'lambdaL1',0);
[nInstances,nVars] = size(X);

if lambdaL1 == 0
    % logisitic L2 regression
    funObj = @(w)LogisticLoss(w,X,y);
    lambda = lambdaL2*ones(nVars, 1);
    lambda(1) = 0;
    options.display= 0;
    w = minFunc(@penalizedL2, zeros(nVars, 1), options, funObj, lambda);
elseif lambdaL2 == 0
    % Least Squares w/ L1 Prior (Lasso)
    w = LassoShooting(X,y,lambdaL1,'verbose',verbose);
end

model.weights = w;
model.predictFunc = @(model,X)sign(X*model.weights);
model.errFunc = @(model,X,y)sum(sign(X*model.weights)~=y)/length(y);
model.lossFunc = @(model,X,y)sum(sign(X*model.weights)~=y)/length(y);
