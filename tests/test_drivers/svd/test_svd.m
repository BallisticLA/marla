%{
    Main routine for running correctness tests of randomized 
    singular value decomposition algorithms (SVD1).

    If all tests pass, no output is present.

    Currently, avoids using logging parameter.

    Important note: 
    Before running tets, use the following:
    
    addpath('../../../utils') - for MatrlaRandStream.m
    addpath('../../../comps/qb') - for different versions of QB algorithm.
    addpath('../../../comps/rangefinders') - for different versions of
    rangefinders.
    addpath('../../../drivers/') - for access to svd1.
%}

% Add all paths here
tsvd1 = TestSVD1();
tsvd1.test_fr();
tsvd1.test_fp_inexact();
tsvd1.test_fp_exact();