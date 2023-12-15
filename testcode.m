addpath('.\func')

%% Create n*n matrix with rank r

n = 10000;
r = 10;
s = 1;

A_0                     = Rank_r_Matrix(n,r,s); % Generate n*n matrix with rank r and about s*100000 nz elements
[A, NZ_col, NZ_row]     = ReducedMatrix(A_0);   % Annihilate all zero rows and columns

figure
subplot(1,2,1)
spy(A_0)
title('Sparsity pattern of A_0')

subplot(1,2,2)
spy(A)
title('Sparsity pattern of A')


%% Compute rref of reduced matrix A
tic
[R_0, jb] = rref(A);
cputime.rref = toc;

%% CR factorization
tic;
[V, D, C, R] = CR(A, R_0, jb);
V = Eigenvector_A_0(A_0, C * V);
cputime.CR = toc;
residual.CR.mtx = mtxdiff(A, C * R);
residual.CR.eig = max(vecnorm(A_0 * V -  V * D));

%% QR factorization
tic;
[V, D, Qr, R] = QR(A);
V = Eigenvector_A_0(A_0, Qr * V);
cputime.QR = toc;
residual.QR.mtx = mtxdiff(A, Qr * R);
residual.QR.eig = max(vecnorm(A_0 * V -  V * D));

%% QR and RREF

tic
fast_R_0 = frref(A);
cputime.Fast_RREF = toc;
residual.Fast_RREF = mtxdiff(R_0,fast_R_0);

%% Krylov method
tic
[V, D] = eigs(A,r);
V = Eigenvector_A_0(A_0, V);
residual.Krylov.eig = max(vecnorm(A_0 * V - V * D));
cputime.Krylov = toc;
