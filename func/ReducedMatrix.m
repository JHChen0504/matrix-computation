function [A ,NZ_cols, NZ_rows] = ReducedMatrix(A)

NZ_cols = max(A ~= 0, [], 1);
NZ_cols = full(NZ_cols);    
A = A(NZ_cols,NZ_cols);

NZ_rows = max(A ~= 0, [], 2);
NZ_rows = full(NZ_rows);  
A = A(NZ_rows,NZ_rows);