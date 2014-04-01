% A simple recommendation system
% Written by Emtiyaz, EPFL, pl8787
% Modified on March 11, 2014

function recommend0(rating_file, movie_file)
    seed = 1;

    % load data
    load(rating_file);
    load(movie_file);
    Y = ratings_train;
    Ytest = ratings_test;
    X = movieMetaData;
    [M,N] = size(Y);

    % Make binary data from Y
    % We will treat ratings with 5 or 1 as binary labels,
    % and treat everything else as missing values
    Y(Y==1) = -1;
    Y(Y==5) = 1;
    Y(Y>1) = 0;
    Ytest(Ytest==1) = -1;
    Ytest(Ytest==5) = 1;
    Ytest(Ytest>1) = 0;

    % Set hyper-params for a GP model
    % There are three (right now set to a suboptimal value)
    logScale = log(1); % GP prior parameters
    logSigma = log(1);
    sig2 = 0.01; % likelihood variance

    % We will use zero mean and covSEiso Kernel
    mu = meanFuncMultiClass({@meanZero}, 2, [], X);
    Sigma = covFuncMultiClass({@covSEiso}, 2, [logScale logSigma], X);

    % Gaussian process regression
    errTr = [];
    errTe = [];
    for n = 1:N
      I = find(Y(:,n));
      Ite = find(Ytest(:,n));
      if size(I,1)==0 || size(Ite,1)==0
          continue
      end
      yn = Y(I,n);
      G = inv(sig2*eye(length(I)) + Sigma(I,I));
      e = G*(yn-mu(I));
      % prediction and errors
      pHat = mu(I) + Sigma(I,I)*e;
      errTr = [errTr; Y(I,n) - pHat];
      pHat = mu(Ite) + Sigma(Ite,I)*e;
      errTe = [errTe; Ytest(Ite,n) - pHat];
    end
    fprintf('Train Error %.4f Test Error %.4f\n', sqrt(mean(errTr.^2)), sqrt(mean(errTe.^2)));
end


