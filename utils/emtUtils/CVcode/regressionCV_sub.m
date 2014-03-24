function [minLambda,minScore] = regressionCV_sub(X,y,options)

if nargin < 3
    options = [];
end

[verbose,trainFunc,lambdaParam,lambdaValues,scaled,prune] = myProcessOptions(options,...
    'verbose',1,'trainFunc',@regressionL2,...
    'CVparam','lambdaL2','CVvalues',2.^[10:-1:-10],...
    'scaled',1,'prune',0);
options.verbose = 0;

minScore = inf;
minLambda = [];
nValues = length(lambdaValues);

if prune < 2

    % Exhaustive Search (prune == 0)
    % or:
    % Exhaustive Search with convex pruning (prune == 1)

    for i = 1:nValues
        % Get next value of lambda
        if prune
            if i == 1
                legal = ones(nValues,1);
                scores = inf*ones(nValues,1);
            end
            candidates = find(legal==1);
            lambdaInd = candidates(ceil(rand*length(candidates)));
        else
            lambdaInd = i;
        end

        [testScore,minScore,minLambda] = CVscore(lambdaValues(lambdaInd),minScore,minLambda,trainFunc,X,y,options,lambdaParam);

        % Prune Parameters
        if prune
            scores(lambdaInd) = testScore;
            computed = find(~isinf(scores));
            [minVal minInd] = min(scores);
            maxLower = max(computed(computed < minInd));
            minUpper = min(computed(computed > minInd));
            legal(lambdaInd) = 0;
            legal(1:maxLower) = 0;
            legal(minUpper:end) = 0;
            if sum(legal==1) == 0
                break;
            end
            lambdaValues(legal==1);
        end
    end
elseif prune == 2
    % Greedy Search (assumes no ties!)
    lambdaInd = floor(nValues/2);

    [testScore,minScore,minLambda] = CVscore(lambdaValues(lambdaInd),minScore,minLambda,trainFunc,X,y,options,lambdaParam);

    % Compute Score with the higher value
    testScoreU = inf;
    if lambdaInd < nValues
        [testScoreU,minScore,minLambda] = CVscore(lambdaValues(lambdaInd+1),minScore,minLambda,trainFunc,X,y,options,lambdaParam);
    end

    % Compute Score with the lower value
    testScoreL = inf;
    if lambdaInd > 1
        [testScoreL,minScore,minLambda] = CVscore(lambdaValues(lambdaInd-1),minScore,minLambda,trainFunc,X,y,options,lambdaParam);
    end

    if testScoreL < testScore && testScoreL < testScoreU
        dir = -1;
        lambdaInd = lambdaInd-1;
        testScoreNew = testScoreL;
    elseif testScoreU < testScore && testScoreU < testScore
        dir = 1;
        lambdaInd = lambdaInd+1;
        testScoreNew = testScoreU;
    else
        lambdaInd = inf;
        testScoreNew = inf;
    end

    while testScoreNew < testScore && lambdaInd >= 2 && lambdaInd <= nValues-1
        testScore = testScoreNew;
        lambdaInd = lambdaInd+dir;
        [testScoreNew,minScore,minLambda] = CVscore(lambdaValues(lambdaInd),minScore,minLambda,trainFunc,X,y,options,lambdaParam);
    end
else
    % 1D Line Search
    lambda1 = lambdaValues(floor(nValues/2));
    lambda3 = lambdaValues(floor(nValues/2)+1);
    lambda2 = (lambda1+lambda3)/2;
    funObj = @(lambda)getScore(lambda,trainFunc,X,y,options,lambdaParam);
    [minLambda,minScore] = brent(funObj,lambda1,lambda2,lambda3,funObj(lambda2));
end

fprintf('minScore = %.3f, min value of %s = %.3f\n',minScore,lambdaParam,minLambda);

end

%%
function [score,minScore,minLambda] = CVscore(lambda,minScore,minLambda,trainFunc,X,y,options,lambdaParam)

score =getScore(lambda,trainFunc,X,y,options,lambdaParam);

% Update max
if score < minScore
    minScore = score;
    minLambda = lambda;
end
end

function [score] = getScore(lambda,trainFunc,X,y,options,lambdaParam)

if lambda < 0 
    score = 1e100;
    return;
end

fprintf('%s = %.3f',lambdaParam,lambda);

nTrain = size(X,1);

score = 0;
for split = 1:2

    % Make Split
    if split == 1
        trainNdx = 1:floor(nTrain/2);
        testNdx = floor(nTrain/2)+1:nTrain;
    else
        trainNdx = floor(nTrain/2)+1:nTrain;
        testNdx = 1:floor(nTrain/2);
    end

    % Set Value
    options=setfield(options,lambdaParam,lambda);

    % Train
    model = trainFunc(X(trainNdx,:),y(trainNdx),options);

    % Compute Test Score
    if isfield(model,'lossFunc')
        score = score + model.lossFunc(model,X(testNdx,:),y(testNdx));
    else
        yhat = model.predictFunc(model,X(testNdx,:)); 
        score = score + sum(abs(y(testNdx) - yhat(:)));
    end
end

if isfield(model,'lossFunc')
fprintf(' (loss = %.3f)\n',score);
else
fprintf(' (CVerr = %.3f)\n',score);
end
end
