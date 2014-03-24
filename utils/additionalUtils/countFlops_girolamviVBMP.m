function nflops = countFlops_girolamviVBMP(N,K,iters)

%K = create_kernel_no_scaling(X,XKernel_Type,Theta,p);
% addflops(countFlop_create_kernel_no_scaling(N,D));
%addflops(countFlops_covSEiso(N,D));
%K = K + eye(N)*SMALL_NOS;
%addflops(N)
% iK = inv(K + In);
nflops = N + flops_inv(N);
% Ki = K*iK;
nflops = nflops + flops_mul(N,N,N);
% for iters: optim bound (for all iters)
nflops = [nflops; zeros(iters-1,1)] + ones(iters,1)*countFlops_vbmpBound(N,K);
return;

function nflops = countFlops_vbmpBound(N,K)

flops(0);
%for k = 1:C, M(:,k) = Ki*Y(:,k); end
addflops(K*flops_mul(N,N,1));
% Loop N: [a,z] =  tmean(M(n,:)',t(n),Nos_Samps_TG);
addflops(N*countFlops_tmean(K))
% lower_bound = lower_bound + safelog(z);
addflops(N)
%iK = inv(K + In);
addflops(N + flops_inv(N));
% Ki = K*iK;
addflops(flops_mul(N,N,N));
% lower_bound = lower_bound +...
% 0.5*N*C -0.5*sum(diag(M'*inv(K)*M))
addflops( 1 + flops_inv(N) + flops_mul(N,N,K) + flops_mul(K,N,K) + K);
% -0.5*C*trace(iK) +0.5*C*logdet(iK);
addflops(flops_det(N) + 1 + N + 2);

nflops = flops;
return;

function nflops = countFlops_tmean(K)
flops(0)
N = 1000;
%u = randn(Nsamps,1);
addflops(flops_randnorm(N));
%t = m(indexMax).*ones(K,1) - m;
addflops(K);
%s = repmat(u,1,K-1) + repmat(t,1,Nsamps)';
addflops(N*K)
% safenormcdf(s')  %thresh + o = 0.5 * erfc(-z ./ sqrt(2));
addflops(K*N + 2*(K-1)*N);
% z = mean(prod(~,1) );
addflops((K-1)*N + K + flops_div);
%loop K: sr = repmat(u,1,K) + repmat(tr,1,Nsamps)';
addflops(K*K*N);
% a1 = u' + m(indexMax)- m(r)
addflops(K*N)
% a2 = safenormcdf(sr')
addflops(K*((K-1)*N + 2*(K-1)*N));
% a3 = safenormpdf(a1) % thresh + o = o.*exp(-((x-m).^2)./(2*(s.^2)));
addflops(K*2*N); % thresh
addflops(K*(2*N + N*flops_div + N*flops_exp));
% a4 = prod(a2,1)
addflops(K*(K-1)*N);
% a5 = mean(a4+ a3)
addflops(K*(N + flops_div));

nflops = flops;
return;

function nflops = countFlop_create_kernel_no_scaling(N,D)

flops(0);
% dist = sum((X.^2)*T,2)*ones(1,ny) +ones(nx,1)*sum((Y.^2)*T,2)' - 2*(X*T*Y');
addflops(N*D + flops_mul(N,D,D) + flops_row_sum(N,D) + flops_mul(N,1,N));
addflops(N*D + flops_mul(N,D,D) + flops_row_sum(N,D) + flops_mul(N,1,N));
addflops(flops_mul_(D,D,N)+ flops_mul(N,D,N));
% K	= exp(-dist(X1,X2,T));
addflops(N^2*flops_exp);

nflops = flops;
return;

function nflops = countFlops_covSEiso(N,D)

flops(0);
%ell = exp(hyp(1)); sf2 = exp(2*hyp(2));
addflops(2*flops_exp)
%K = sq_dist(x'/ell,z'/ell);
addflops(N*D + flops_mul(N,D,N));
%K = sf2*exp(-K/2);
addflops(N^2 + N^2*flops_exp);

nflops = flops;
return;
