%{
    Main routine for running correctness tests of randomized 
    QB algorithms (QB1, QB2, QB3, QB4).

    If all tests pass, no output is present.

    Currently, avoids using logging parameter.

    Important note: 
    Before running tets, use the following:
    
    addpath('../../utils') - for MatrlaRandStream.m
    addpath('../../comps/qb') - for different versions of QB algorithm.
    addpath('../../comps/rangefinders') - for different versions of
    rangefinders.
    addpath('../matrix_generators/') - for matrix generators
%}
% QB-backed algorithm
addpath('../../utils') 
addpath('../../comps/qb') 
addpath('../../comps/rangefinders') 
addpath('../matrix_generators/')

tqb = TestQB1();
tqb.run_var_rank();