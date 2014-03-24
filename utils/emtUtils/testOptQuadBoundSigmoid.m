% test optComputeBoundSigmoid code
clear all
yRange = [-5:.1:5];
mean_ = 1;
variance = 10;
% compute exact value
logSig = log(1+exp(yRange));
% compute lower bound
[xi, A, b, c] = optQuadBoundSigmoid(mean_, variance);
lb = 0.5*yRange.^2*A - b*yRange + c;

% plot
plot(yRange, logSig,'linewidth',2);
hold on
plot(yRange, lb, 'r','linewidth',2);
legend('log(1+exp(z))','quadBound');
title(sprintf('Lower bound for mean = %d, variance = %1.2f',mean_, variance));
xlabel('z');
