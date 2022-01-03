classdef TestSVD1 < TestSVDecomposer
%{
    Test class, with functions that represent runs of an rsvd slgorithm
    with varying input parameters. Uses TestSVDecomposer as a
    subclass to run unit tests in batches. 
%}  
    methods
        function obj = TestSVD1()
            obj = obj@TestSVDecomposer();
        end
        
        % Fixed rank tets - target tolerance is set to NaN, so the
        % algorithm terminates upon reaching rank k.
        function [] = test_fr(obj)
            num_passes = 3;
            block_size = 5;
            logging.depth = 0;
            logging.span = 0;
            alg = @(A_, k_, tol_, over_, seed_, logging) svd1(A_, k_, tol_, over_,...
                num_passes, block_size, seed_, logging);

            rank = 15;
            seed = 99;
            alg_tol = 1e-8;
            % Tall matrix
            [U, s, Vt] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            sth = SvdTestHelper(U * diag(s) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank - 10, NaN, alg_tol, 0, logging);
            obj.run_batch(sth, alg, rank - 10, NaN, alg_tol, 5, logging);
            obj.run_batch(sth, alg, rank - 3, NaN, alg_tol, 5, logging);

            % Wide matrix
            [U, s, Vt] = gen_test_mat_factors(50, 200, rank, NaN, seed);
            sth = SvdTestHelper(U * diag(s) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank - 10, NaN, alg_tol, 0, logging);
            obj.run_batch(sth, alg, rank - 10, NaN, alg_tol, 5, logging);
            obj.run_batch(sth, alg, rank - 3, NaN, alg_tol, 5, logging);
        end
        
        % Test with a ~bad quality approximation - algorithm  tolerance is 
        % set to a relatively high value.
        function [] = test_fp_inexact(obj)
            num_passes = 2;
            block_size = 5;
            logging.depth = 0;
            logging.span = 0;            
            alg = @(A_, k_, tol_, over_, seed_, logging) svd1(A_, k_, tol_, over_,...
                num_passes, block_size, seed_, logging);
            seed = 11111;
            rel_err = 0.001;
            rank = 50;

            [U, s, Vt] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            sth = SvdTestHelper(U * diag(s) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0, logging);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2, logging);

            [U, s, ~] = gen_test_mat_factors(50, 200, rank, NaN, seed);
            sth = SvdTestHelper(U * diag(s) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0, logging);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2, logging);
        end

        % Test with algorithm tolerance lower tan test tolrance - good
        % approximation.
        function [] = test_fp_exact(obj)
            num_passes = 2;
            block_size = 5;
            logging.depth = 0;
            logging.span = 0;            
            alg = @(A_, k_, tol_, over_, seed_, logging) svd1(A_, k_, tol_, over_,...
                num_passes, block_size, seed_, logging);
            seed = 1512;
            rank = 15;
            rel_err = 1e-12;

            [U, s, Vt] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            sth = SvdTestHelper(U * diag(s) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0, logging);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2, logging);

            [U, s, Vt] = gen_test_mat_factors(50, 200, rank, NaN, seed);
            sth = SvdTestHelper(U * diag(s) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0, logging);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2, logging);

            [U, s, Vt] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            sth = SvdTestHelper(U * diag(s) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0, logging);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2, logging);
        end
    end
end