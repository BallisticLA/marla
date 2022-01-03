function [A_bad] = gen_simp_mat(n_rows, n_cols, scale, s)
%{
    QR - based routine for random matrix generation.

    Input
    -----
    n_rows, n_cols : int
        Matrix dimensions

    scale :

    s : int or randstream
        Random seed.
 
     Output
     ------
     A_bad : matrix
        Size n_rwos by n_cols

    Important note: 
    ---------------
    Before calling this routine, use:
    addpath('../../../utils/');
%}
    s = MarlaRandStream(s);
    A = randn(s, n_rows, n_cols);
    [QA, RA] = qr(A);
    damp = diag(arrayfun(@(x) inv(x), sqrt(1 + scale * (1:n_cols))));
    RA = RA * damp';
    A_bad = QA * RA;
end