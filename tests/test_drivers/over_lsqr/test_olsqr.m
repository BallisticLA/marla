%{
    Main routine for running correctness tests of randomized 
    overdetermined least squares algorithms (SPO1 and SPO2).

    If all tests pass, no output is present.

    Currently, avoids using logging parameter.

    Important note: 
    Before running tets, use the following:
    
    addpath('../../../utils') - for MatrlaRandStream.m
    addpath('../../../comps/qb') - for different versions of QB algorithm.
    addpath('../../../comps/rangefinders') - for different versions of
    rangefinders.
    addpath('../../../drivers/least_squares') - for different versions of evd.
%}

addpath('../../../utils/sketching_operators/');
addpath('../../../utils');
addpath('../../matrix_generators/');
addpath('../../../comps/determiter/')
addpath('../../../drivers/least_squares');

spo1 = TestSPO1();
spo1.consistent_tall();
spo1.consistent_lowrank();
spo1.consistent_square();
spo1.inconsistent_orthog();
spo1.inconsistent_gen();

spo3 = TestSPO3();
spo3.consistent_tall();
spo3.consistent_square();
spo3.inconsistent_orthog();
spo3.inconsistent_gen();

