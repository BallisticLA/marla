function [U, spectrum, Vt] = gen_test_mat_factors(m, n, k, spectrum, s)
%{
    A routine for generating a random matrix, outputting it in a form of a
    Singular Value Dcomposition

    Input
    -----
    m, n : int
        Matrix dimensions

    k : int
        Matrix rank

    spectrum : scalar or vector
        Used for generating the spectrum of a matrix.

    s : int or randstream
        Random seed.
 
     Output
     ------
     U : matrix
         Orthonormal
 
     S : vetcor
         Size 1 by k vetcor of singular values

     V : matrix
         Orthonormal

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
    Vt = V';
    if isscalar(spectrum)
        spectrum = abs(randn(s, 1, k));
        spectrum = sort(spectrum,'descend');
    end
end