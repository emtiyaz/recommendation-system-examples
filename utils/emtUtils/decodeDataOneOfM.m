function dataNew = decodeDataOneOfM(data, nClass, classSize);
% returns one of M encoding of the data according to number of classes nClass.
% Data should be #dim x #measurments
% if classSize is 'M+1' then only first nClass-1 bits are returned.
% if classSize is 'M' then all the bits are returned
  if nargin==2
    classSize = 'M+1';
  end

  dataNew = [];
  for d = 1:length(nClass)
    if strcmp(classSize,'M+1')
      % collect the dth data
      M = nClass -1;
      idx = sum(M(1:d-1))+1:sum(M(1:d));
      y = data(idx,:); 
      lastRow = sum(y,1);
      obs = ~isnan(lastRow);
      lastRow(obs) = ~lastRow(obs) ;
      y = [y; lastRow];
      % convert to integers
      dataNew(d,:) = sum(bsxfun(@times, y, [1:nClass(d)]'),1);
    else
      % collect the dth data
      M = nClass;
      idx = sum(M(1:d-1))+1:sum(M(1:d));
      y = data(idx,:); 
      % convert to integers
      dataNew(d,:) = sum(bsxfun(@times, y, [1:nClass(d)]'),1);
    end
  end

