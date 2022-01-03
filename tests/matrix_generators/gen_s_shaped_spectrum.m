function [A, S] = gen_s_shaped_spectrum(m, n, k, s)
%{
    A routine for random matrix generation, additionally outputting its
    spectrum in a form of a diagonal matrix.
    
    Singular values have an S-shaped decay when plotted.

    Input
    -----
    m, n : int
        Matrix dimensions

    k : int
        Matrix rank

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
    spectrum = 0.0001+1./(1 + exp((1:k)-30));
    [A, S] = gen_test_mat(m, n, k, spectrum, s);
end