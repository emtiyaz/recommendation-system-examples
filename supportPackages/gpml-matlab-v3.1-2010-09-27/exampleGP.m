clear all;
y = 2*double(rand(10,1)>0.5)-1;
X = rand(10,15);
ys = 2*double(rand(5,1)>0.5)-1;
Xs = rand(5,15);
meanfunc = @meanZero; 
covfunc = @covSEiso;   hyp.cov = log([0.5 0.5]);
likfunc = @likLogistic;
[ymu ys2 fmu fs2 lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, X, y, Xs, ys);
