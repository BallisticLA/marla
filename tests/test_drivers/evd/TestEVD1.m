classdef TestEVD1 < TestEVDecomposer
%{
    Test class for EVD1. Uses TestEVDecomposer as a subclass to run unit
    tests in batches.
%}  
    
    methods
        function obj = TestEVD1(psd)
            obj = obj@TestEVDecomposer(psd);
        end
        
        % Fixed rank tets - target tolerance is set to NaN, so the
        % algorithm terminates upon reaching rank k.
        function [] = test_fr(obj)
            addpath('../../matrix_generators/');
            addpath('../../../drivers/evd/');
            num_passes = 3;
            block_size = 2;
            logging.depth = 0;
            logging.span = 0;
            alg = @(A_, k_, tol_, over_, seed_, logging) evd1(A_, k_, tol_, over_,...
                num_passes, block_size, seed_, logging);

            rank = 15;
            seed = 7;
            [U, s, ~] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            eth = EigTestHelper.convert(U, s, obj.PSD, seed);

            obj.run_batch(eth, alg, rank - 10, NaN, 1e-8, 0, logging);
            obj.run_batch(eth, alg, rank - 10, NaN, 1e-8, 5, logging);
            obj.run_batch(eth, alg, rank - 3, NaN, 1e-8, 5, logging);
        end

        % Test with a ~bad quality approximation - algorithm  tolerance is 
        % set to a relatively high value.
        function [] = test_fp_inexact(obj)
            addpath('../../matrix_generators/');
            addpath('../../../drivers/evd/');
            num_passes = 2;
            block_size = 2;
            logging.depth = 0;
            logging.span = 0;
            alg = @(A_, k_, tol_, over_, seed_, logging) evd1(A_, k_, tol_,...
                over_, num_passes, block_size, seed_, logging);
            seed = 0;

            rank = 15;
            [U, s, ~] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            eth = EigTestHelper.convert(U, s, obj.PSD, seed);

            % TODO: update tests so we can verify returned matrix
            % has a reasonably low rank.
            rel_err = 0.25;
            n = size(U, 1);
            obj.run_batch(eth, alg, n, rel_err, 1e-8, 0, logging);
            obj.run_batch(eth, alg, n, rel_err, 1e-8, 3, logging);
        end

        % Test with algorithm tolerance lower tan test tolrance - good
        % approximation.
        function [] = test_fp_exact(obj)
            num_passes = 2;
            block_size = 2;
            logging.depth = 0;
            logging.span = 0;
            alg = @(A_, k_, tol_, over_, seed_, logging) evd1(A_, k_, tol_,...
                over_, num_passes, block_size, seed_, logging);
            seed = 0;

            rank = 15;
            [U, s, ~] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            eth = EigTestHelper.convert(U, s, obj.PSD, seed);

            rel_err = 1e-12;
            obj.run_batch(eth, alg, rank, rel_err, 1e-8, 0, logging);
            obj.run_batch(eth, alg, rank, rel_err, 1e-8, 3, logging);
        end
    end
end

