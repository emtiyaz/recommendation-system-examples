% last two options don't work very well right now
nInstances = 500;
nVars = 100;
X = [ones(nInstances,1) randn(nInstances,nVars-1)];
w = randn(nVars,1);
w(nVars/2+1:end) = 0;
y = sign(X*w);

fprintf('Cross-validation by enumeration...\n');
options = struct('trainFunc',@classificationL2,'CVparam','lambdaL2','CVvalues',2.^[10:-1:-10]);
model = regressionCV(X,y,options);
break

pause;

fprintf('Cross-validation by enumeration (with pruning)...\n');
options.prune = 1;
model = regressionCV(X,y,options);
pause;

fprintf('Cross-validation by greedy...\n');
options.prune = 2;
model = regressionCV(X,y,options);
pause;

% Cross-validation with line-minimization
fprintf('Cross-validation by line minimization...\n');
options.prune = 3;
options.CVvalues = [0 1000];
options = struct('trainFunc',@classificationL2,'CVparam','lambdaL2','CVvalues',[0 1000],'prune',3);
model = regressionCV(X,y,options);
