function val = logMvnPdf(y, mu, precMat)

  D = size(mu,2);
  for i = 1:size(y,1)
    diff = y(i,:)-mu;
    val(i) = -0.5*D*log(2*pi) + 0.5*logdet(precMat) - 0.5*diff*precMat*diff';
  end
