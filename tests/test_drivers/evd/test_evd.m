%{
    Main routine for running correctness tests of randomized 
    Eigendecomposition algorithms (EVD1 and EVD2).

    If all tests pass, no output is present.

    Currently, avoids using logging parameter.

    Important note: 
    Before running tets, use the following:
    
    addpath('../../../utils') - for MatrlaRandStream.m
    addpath('../../../comps/qb') - for different versions of QB algorithm.
    addpath('../../../comps/rangefinders') - for different versions of
    rangefinders.
    addpath('../../../drivers/evd') - for different versions of evd.
%}
% QB-backed algorithm
tevd1 = TestEVD1(false);
tevd1.test_fr();
tevd1.test_fp_inexact();
tevd1.test_fp_exact();

% Nystrom, for PSD matrices.
tevd2 = TestEVD2();
tevd2.test_fr();