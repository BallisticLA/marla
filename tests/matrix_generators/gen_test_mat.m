function [A, S] = gen_test_mat(m, n, k, spectrum, s)
%{
    A routine for random matrix generation, additionally outputting its
    spectrum in a form of a diagonal matrix.

    Input
    -----
    m, n : int
        Matrix dimensions

    k : int
        Matrix rank

    spectrum : int or vector
        Used for generating the spectrum of a matrix.

    s : int or randstream
        Random seed.
 
     Output
     ------
     A : matrix
        Size m by n random matrix of rank k having a spctrum based on
        provided specification.

     S : matrix
        Diagonal matrix of singular values of A.

    Important note: 
    ---------------
    Before calling this routine, use:
    addpath('../../../utils/');
%}
    s = MarlaRandStream(s);
    Buf = randn(s, m, k);
    [U, ~] = qr(Buf, 0);
    Buf = randn(s, n, k);
    [V, ~] = qr(Buf, 0);
    if isscalar(spectrum)
        spectrum = abs(randn(s, 1, k));
        spectrum = sort(spectrum,'descend');
    end
    S = spdiags(spectrum', 0, k, k);
    A = U * S * V';
end