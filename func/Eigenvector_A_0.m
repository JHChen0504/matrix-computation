function V = Eigenvector_A_0(A_0,V)

[~, NZ_col, NZ_row] = ReducedMatrix(A_0);

A1 = A_0(:,NZ_col);

VV = sparse(length(NZ_row),size(V,2));

VV(NZ_row,:) = V;

V = A1*VV;