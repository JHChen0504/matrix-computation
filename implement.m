clear;clc;close all;
rng(1)

% Create n*n matrix with rank r

n = 1000;
r = 5;
P = randperm(n);

A = [speye(r) sparse(round(normrnd(0,0.2,r,n-r)))];
A = A(:,P);
X = sparse(n,n);

for i = 1:r
    X = X + A(i,:).'*A(i,:);
end

figure
spy(X)

%% CR factorization

cr = tic;

[R0,jb] = rref(X);
rank_X = length(jb);

C = X(:,jb);
R = R0(1:rank_X,:);
Xrc = R*C;
cpu.cr = toc(cr);

%% CR factorization with fast rref

fcr = tic;

[R0,jb] = frref(X);
rank_X = length(jb);

Cf = X(:,jb);
Rf = R0(1:rank_X,:);
Xrcf = Rf*Cf;
cpu.fcr = toc(fcr);

%% SVD 

svu = tic;

X = full(X);
[U,S,V] = svd(X);
Uk = U(:,1:rank_X);
Sk = diag(S);
Sk = diag(Sk(1:rank_X));
Vk = V(:,1:rank_X);

Xsvu = Sk*Vk'*Uk;
cpu.svu = toc(svu);

%% compute eigenvalue

tx = tic;
[V , D] = eig(X);
cpu.tx = toc(tx);

tcr = tic;
[V_cr,D_cr] = eig(Xrc);
cpu.tcr = toc(tcr);

Xrcf = full(Xrcf);
ftcr = tic;
[V_fcr,D_fcr] = eig(Xrcf);
cpu.ftcr = toc(ftcr);

tsvd = tic;
[V_svd,D_svd] = eig(Xsvu);
cpu.tsvd = toc(tsvd);

%% Computationl Error

fprintf('Computationl Error : \n\n')

fprintf('||X-CR||      = %2.4e\n',norm(X-C*R))
fprintf('||X-CfRf||    = %2.4e\n',norm(X-Cf*Rf))
fprintf('||X-USV''||    = %2.4e\n',norm(X-U*S*V'))
fprintf('||X-UkSkVk''|| = %2.4e\n\n%',norm(X-Uk*Sk*Vk'))

fprintf('Computationl Error eig : \n\n')

eX = maxk(diag(D),rank_X,'ComparisonMethod','abs');

fprintf('eigenvalue(X-CR)      = %2.4e\n',sqrt(sum((eX-maxk(abs(diag(D_cr)),rank_X)).^2)))
fprintf('eigenvalue(X-CfRf)    = %2.4e\n',sqrt(sum((eX-maxk(abs(diag(D_fcr)),rank_X)).^2)))
fprintf('eigenvalue(X-UkSkVk'') = %2.4e\n\n%',sqrt(sum((eX-maxk(abs(diag(D_svd)),rank_X)).^2)))


%% Elapsed CPU Time for eig

fprintf('Elapsed CPU Time computing eig : \n\n')

fprintf('eig(X)            = %fs\n',cpu.tx)
fprintf('eig(RC)           = %fs\n',cpu.tcr)
fprintf('eig(RfCf)         = %fs\n',cpu.ftcr)
fprintf('eig(SkVk''Uk)      = %fs\n\n',cpu.tsvd)

%% Elapsed CPU Time for whole process

fprintf('Elapsed CPU Time for whole process : \n\n')

fprintf('eig(X)            = %fs\n',cpu.tx)
fprintf('eig(RC)           = %fs\n',cpu.cr+cpu.tcr)
fprintf('eig(RfCf)         = %fs\n',cpu.fcr+cpu.ftcr)
fprintf('eig(SkVk''Uk)      = %fs\n\n',cpu.svu+cpu.tsvd)


