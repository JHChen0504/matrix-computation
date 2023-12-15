function X = Rank_r_Matrix(n,r,sparsity)

arguments
    n double {mustBeInteger,mustBePositive,mustBeFinite}
    r double {mustBeInteger,mustBePositive,mustBeFinite}
    sparsity double = 1
end

if r>n
    error('Invalid input rank. Rank value must be less than n = %d.',n)
end

s = sqrt(sparsity/ (0.1*r))*(0.1^log10(n)) *100;

P = randperm(n);

A = [speye(r) sprand(r , n - r, s)];
A = A(:,P);

C = sprand(n - r , r , s);

C(C ~= 0) = 2*C(C ~= 0)-1;

X = [A ; C * A];
X = X(P,:);