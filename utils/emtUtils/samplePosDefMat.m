function C = samplePosDefMat(L);
% sample a posdef matrix with condition number 10
% Written by Emtiys, CS, UBC

    C = randn(L,L);
    [Q,R] = qr(C);
    C = Q*diag(linspace(0.1, 1, L))*Q';

