function [trainData, testData, idx] = splitData(yc, yd, ratio) 
% splits data into training and testing set
% yc is the continuous data, yd is discrete data
% ratio is the split ratio

  [Dc,Nc] = size(yc);
  [Dd,Nd] = size(yd);
  N = max(Nc,Nd);
  nTrain = ceil(ratio*N);
  idx = randperm(N);
  if Dc>0
    testData.continuous = yc(:,idx(nTrain+1:end));
    trainData.continuous = yc(:,idx(1:nTrain));
  else
    testData.continuous = [];
    trainData.continuous = [];
  end
  if Dd>0
    testData.discrete = yd(:,idx(nTrain+1:end));
    trainData.discrete = yd(:,idx(1:nTrain));
  else
    testData.discrete = [];
    trainData.discrete = [];
  end
 
