function [model] = regressionCV(X,y,options)
% w = regressionCV(X,y)
%
% Runs sub-classifier using cross-validation to select parameters
%
% Scaled should be 1 for parameters that depend on the number of data
% points (like regularization), and 0 otherwise
%%

if nargin < 3
    options = [];
end

[verbose,trainFunc,lambdaParam,lambdaValues,scaled,prune] = myProcessOptions(options,...
    'verbose',1,'trainFunc',@regressionL2,...
    'CVparam','lambdaL2','CVvalues',2.^[10:-1:-10],...
    'scaled',1,'prune',0);

%% Determine Best Parameters by 1D Search
[minLambda,minScore] = regressionCV_sub(X,y,options);


%% Now train final model
if scaled
    minLambda = minLambda*2;
end
fprintf('Setting cross-validated parameter to %f\n',minLambda);
options=setfield(options,lambdaParam,minLambda);
model = trainFunc(X,y,options);

end


