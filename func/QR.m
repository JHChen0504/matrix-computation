function [V,D,Qr,R] = QR(A) 
    [Q, R_0] = qr(A);
    NZ_rows = max(R_0 ~= 0, [], 2);
    NZ_rows = full(NZ_rows); % probably more efficient
    NZ_rows(~NZ_rows) = [];
    rank_QR = length(NZ_rows);
    Qr = Q(:,1:rank_QR);
    R = R_0(NZ_rows,:);
    X = R * Qr;
    [V, D] = eigs(X,rank_QR);
end