function[M, U, S, V] = svd_right_precond(A_ske)
    [U, S, V] = svd(A_ske, 'econ', 'vector');
    U = U(:, 1:rank);
    rank = nnz(S > S(1) * num_cols * eps('double'));
    M = V(:, 1:rank) ./ (S(1:rank)');
end