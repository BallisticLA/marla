function[M, U, S, V] = svd_right_precond(A_ske)
    [U, S, V] = svd(A_ske, 'econ');
    S = diag(S);
    rk = nnz(S > S(1) * size(A_ske, 2) * eps('double'));
    U = U(:, 1:rk);
    M = V(:, 1:rk) ./ (S(1:rk)');
end