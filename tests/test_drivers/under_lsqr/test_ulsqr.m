%{
    Main routine for running correctness tests of randomized 
    underdetermined least squares algorithms (SPU1).

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
addpath('../../../comps/preconditioning/')
addpath('../../../drivers/least_squares');

spu1 = TestSPU1();
spu1.linspace_spec();
spu1.logspace_spec();
spu1.lowrank_linspace_spec();
