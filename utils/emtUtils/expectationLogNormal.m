function out = expectationLogNormal(meanP, precMatP, meanQ, covMatQ)
% compute E(log Normal(meanP, inv(precMatP))) under N(meanQ, covMatQ)
% written by Emtiyaz, CS, UBC
% modified on Feb 09, 2010

  D = length(meanP);
  diff = meanP - meanQ;
  out = -0.5*D*log(2*pi) + 0.5*logdet(precMatP) - 0.5*diff'*precMatP*diff...
        -0.5*trace(precMatP*covMatQ);

