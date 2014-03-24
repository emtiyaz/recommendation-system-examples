function mu = my_softmax(eta)

mu = exp(eta)./sum(exp(eta));
