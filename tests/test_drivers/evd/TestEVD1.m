classdef TestEVD1 < TestEVDecomposer
    
    methods
        function obj = TestEVD1(psd)
            obj = obj@TestEVDecomposer(psd);
        end
        
        function [] = test_fr(obj)
            num_passes = 3;
            block_size = 2;
            alg = @(A_, k_, tol_, over_) evd1(A_, k_, tol_, over_,...
                num_passes, block_size);
            state = rng(0);

            rank = 15;
            [U, s, ~] = gen_test_mat_factors(200, 50, rank, NaN);
            [eth, state] = EigTestHelper.convert(U, s, obj.PSD, state);
            rng(state);

            obj.run_batch(eth, alg, rank - 10, NaN, 1e-8, 0);
            obj.run_batch(eth, alg, rank - 10, NaN, 1e-8, 5);
            obj.run_batch(eth, alg, rank - 3, NaN, 1e-8, 5);
        end

        function [] = test_fp_inexact(obj)
            num_passes = 2;
            block_size = 2;
            alg = @(A_, k_, tol_, over_) evd1(A_, k_, tol_, over_,...
                num_passes, block_size);
            state = rng(0);

            rank = 15;
            [U, s, ~] = gen_test_mat_factors(200, 50, rank, NaN);
            [eth, state] = EigTestHelper.convert(U, s, obj.PSD, state);
            rng(state);
            % TODO: update tests so we can verify returned matrix
            % has a reasonably low rank.
            rel_err = 0.25;
            n = size(U, 1);
            obj.run_batch(eth, alg, n, rel_err, 1e-8, 0);
            obj.run_batch(eth, alg, n, rel_err, 1e-8, 3);
        end

        function [] = test_fp_exact(obj)
            num_passes = 2;
            block_size = 2;
            alg = @(A_, k_, tol_, over_) evd1(A_, k_, tol_, over_,...
                num_passes, block_size);
            state = rng(0);

            rank = 15;
            [U, s, ~] = gen_test_mat_factors(200, 50, rank, NaN);
            [eth, state] = EigTestHelper.convert(U, s, obj.PSD, state);
            rng(state);

            rel_err = 1e-12;
            obj.run_batch(eth, alg, rank, rel_err, 1e-8, 0);
            obj.run_batch(eth, alg, rank, rel_err, 1e-8, 3);
        end
    end
end

