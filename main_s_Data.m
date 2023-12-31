addpath('.\func')

n = 10000;
r = 500;
s = linspace(1,2.5,10);

N = length(s);
cputime.QR = zeros(1,N);
cputime.Krylov = zeros(1,N);
residual.QR.mtx = zeros(1,N);
residual.QR.eig = zeros(1,N);
residual.Krylov.eig = zeros(1,N);

tstart = tic;

for i = 1:N
    A_0 = Rank_r_Matrix(n,r,s(i));
    [A, NZ_col, NZ_row]     = ReducedMatrix(A_0);
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

l = 'Matrix sparsity';
t = ['dimension n = ' num2str(n) ' rank r = ' num2str(r)];

time = figure;
plot(s,cputime.QR,'b-o',s,cputime.Krylov,'r-o','LineWidth',2,'MarkerSize',3)
xlabel(l)
ylabel('CPU time')
legend('QR','Krylov','Location','NW','Fontsize',14)
title(t,'FontSize',14)

Reseig = figure;
plot(s,residual.QR.eig,'b-o',s,residual.Krylov.eig,'r-o','LineWidth',2,'MarkerSize',3)
xlabel(l)
ylabel('Residual of eigenvector')
legend('QR','Krylov','Location','NW','Fontsize',14)
title(t,'FontSize',14)

Resmtx = figure;
plot(s,residual.QR.mtx,'b-o','LineWidth',2,'MarkerSize',3)
xlabel(l)
ylabel('Residual of matrix')
legend('QR','Location','NW','Fontsize',14)
title(t,'FontSize',14)

saveas(time,'s_data_time','png')
saveas(Reseig,'s_data_eig','png')
saveas(Resmtx,'s_data_mtx','png')

