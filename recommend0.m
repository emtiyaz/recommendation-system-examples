% A simple recommendation system
% Written by Emtiyaz, EPFL
% Modified on March 11, 2014

%clear all
seed = 1;

% load data
load('movielens100k.mat');
Y = ratings;
X = movieMetaData;
[M,N] = size(Y);

% Make binary data from Y
% We will treat ratings with 5 or 1 as binary labels,
% and treat everything else as missing values
Y(Y==1) = -1;
Y(Y==5) = 1;
Y(Y>1) = 0;

% Make train and test set
setSeed(seed);
[D N] = size(Y);
num = sum(Y~=0,2);
dd = []; nn = []; yy = [];
for n = 1:N
    On = find(Y(:,n)~=0);
    if length(On)>10 % if user has more than 10 ratings
        ind = unidrnd(length(On));
        d = On(ind);
        if num(d)>10 % movie has more than 10 ratings
            dd = [dd; d];
            nn = [nn; n];
            yy = [yy; Y(d,n)];
        end
    end
end
Ytest = sparse(dd,nn,yy,D,N);
Y(sub2ind([D N], dd, nn)) = 0;

% Set hyper-params for a GP model
% There are three (right now set to a suboptimal value)
logScale = 0; % GP prior parameters
logSigma = 0;
sig2 = .01; % likelihood variance

% We will use zero mean and covSEiso Kernel
mu = meanFuncMultiClass({@meanZero}, 2, [], X);
Sigma = covFuncMultiClass({@covSEiso}, 2, [logScale logSigma], X);

% Gaussian process regression
errTr = [];
errTe = [];
for n = 1:N
  I = find(Y(:,n));
  yn = Y(I,n);
  G = inv(sig2*eye(length(I)) + Sigma(I,I));
  e = G*(yn - mu(I));
  % prediction and errors
  pHat = mu(I) + Sigma(I,I)*e;
  errTr = [errTe; Y(I,n) - pHat];
  Ite = find(Ytest(:,n));
  pHat = mu(Ite) + Sigma(Ite,I)*e;
  errTe = [errTe; Ytest(Ite,n) - pHat];
end
fprintf('Train Error %.4f Test Error %.4f\n', sqrt(mean(errTr.^2)), sqrt(mean(errTe.^2)));


