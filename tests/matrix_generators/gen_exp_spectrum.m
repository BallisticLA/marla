function [A, S] = gen_exp_spectrum(m, n, k, t, s)
%{
    A routine for random matrix generation, additionally outputting its
    spectrum in a form of a diagonal matrix.

    Input
    -----
    m, n : int
        Matrix dimensions

    k : int
        Matrix rank

    t : int
        Controls the dcay of sungular values. The higher the parameter, the
        "slower" the decay is (singhular values ar close together).

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
    spectrum = exp((1 : k) / -t);
    [A, S] = gen_test_mat(m, n, k, spectrum, s);
end