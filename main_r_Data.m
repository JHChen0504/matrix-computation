addpath('.\func')

n = 10000;
r = 500:50:1000;
s = 1;

N = length(r);
cputime.QR = zeros(1,N);
cputime.Krylov = zeros(1,N);
residual.QR.mtx = zeros(1,N);
residual.QR.eig = zeros(1,N);
residual.Krylov.eig = zeros(1,N);

tstart = tic;

for i = 1:N
    A_0 = Rank_r_Matrix(n,r(i),s);
    [A, NZ_col, NZ_row]     = ReducedMatrix(A_0);
    tic;
    [V, D, Qr, R] = QR(A);
    V = Eigenvector_A_0(A_0, Qr * V);
    cputime.QR(i) = toc;
    residual.QR.mtx(i) = mtxdiff(A, Qr * R);
    residual.QR.eig(i) = max(vecnorm(A_0 * V -  V * D));
    tic
    [V, D] = eigs(A,r(i));
    V = Eigenvector_A_0(A_0, V);
    residual.Krylov.eig(i) = max(vecnorm(A_0 * V - V * D));
    cputime.Krylov(i) = toc;
end
toc(tstart)

%% 

l = 'Matrix rank';
t = ['dimension n = ' num2str(n) ' nz = ' num2str(s*100000)];

time = figure;
plot(r,cputime.QR,'b-o',r,cputime.Krylov,'r-o','LineWidth',2,'MarkerSize',3)
xlabel(l)
ylabel('CPU time')
legend('QR','Krylov','Location','NW','Fontsize',14)
title(t,'FontSize',14)

Reseig = figure;
plot(r,residual.QR.eig,'b-o',r,residual.Krylov.eig,'r-o','LineWidth',2,'MarkerSize',3)
xlabel(l)
ylabel('Residual of eigenvector')
legend('QR','Krylov','Location','NW','Fontsize',14)
title(t,'FontSize',14)

Resmtx = figure;
plot(r,residual.QR.mtx,'b-o','LineWidth',2,'MarkerSize',3)
xlabel(l)
ylabel('Residual of matrix')
legend('QR','Location','NW','Fontsize',14)
title(t,'FontSize',14)

saveas(time,'r_data_time','pdf')
saveas(Reseig,'r_data_eig','pdf')
saveas(Resmtx,'r_data_mtx','pdf')

