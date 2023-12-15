function x = mtxdiff(A,B)
    if size(A) ~= size(B)
        error('Input matrix should have same size')
    end
    x = max(max(abs(A-B)));
end