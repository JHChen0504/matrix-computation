function [V,D,C,R] = CR(A,R_0,jb)
    switch nargin
        case 1
            [R_0,jb] = rref(A);
        case 2
            [Non_Zero_Rows, jb] = max(R_0 ~= 0, [], 2);
            Non_Zero_Rows = full(Non_Zero_Rows); % probably more efficient
            jb = jb(Non_Zero_Rows);
    end
    rank_CR = length(jb);
    R = R_0(1:rank_CR,:);
    C = A(:,jb);
    X = R * C;
    X = full(X);
    [V, D] = eigs(X,rank_CR);
end