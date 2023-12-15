addpath('.\func')

n = 5000:500:10000;
r = 500;
s = 1;

N = length(n);
cputime.QR = zeros(1,N);
cputime.Krylov = zeros(1,N);
residual.QR.mtx = zeros(1,N);
residual.QR.eig = zeros(1,N);
residual.Krylov.eig = zeros(1,N);

tstart = tic;

for i = 1:N
    A_0 = Rank_r_Matrix(n(i),r,s);
    [A, NZ_col, NZ_row] = ReducedMatrix(A_0);
    tic;
    [V, D, Qr, R] = QR(A);
    V = Eigenvector_A_0(A_0, Qr * V);
    cputime.QR(i) = toc;
    residual.QR.mtx(i) = mtxdiff(A, Qr * R);
    residual.QR.eig(i) = max(vecnorm(A_0 * V -  V * D));
    tic
    [V, D] = eigs(A,r);
    V = Eigenvector_A_0(A_0, V);
    residual.Krylov.eig(i) = max(vecnorm(A_0 * V - V * D));
    cputime.Krylov(i) = toc;
end
toc(tstart)

%% 

l = 'Matrix dimesion';
t = ['Rank r = ' num2str(r) ' nz = ' num2str(s*100000)];

time = figure;
plot(n,cputime.QR,'b-o',n,cputime.Krylov,'r-o','LineWidth',2,'MarkerSize',3)
xlabel(l)
ylabel('CPU time')
legend('QR','Krylov','Location','NW','Fontsize',14)
title(t,'FontSize',14)

Reseig = figure;
plot(n,residual.QR.eig,'b-o',n,residual.Krylov.eig,'r-o','LineWidth',2,'MarkerSize',3)
xlabel(l)
ylabel('Residual of eigenvector')
legend('QR','Krylov','Location','NW','Fontsize',14)
title(t,'FontSize',14)

Resmtx = figure;
plot(n,residual.QR.mtx,'b-o','LineWidth',2,'MarkerSize',3)
xlabel(l)
ylabel('Residual of matrix')
legend('QR','Location','NW','Fontsize',14)
title(t,'FontSize',14)

saveas(time,'n_data_time','png')
saveas(Reseig,'n_data_eig','png')
saveas(Resmtx,'n_data_mtx','png')
