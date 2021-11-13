function [A_bad] = gen_simp_mat(n_rows, n_cols, scale, s)
    s = MarlaRandStream(s);
    A = randn(s, n_rows, n_cols);
    [QA, RA] = qr(A);
    damp = diag(arrayfun(@(x) inv(x), sqrt(1 + scale * (1:n_cols))));
    RA = RA * damp';
    A_bad = QA * RA;
end