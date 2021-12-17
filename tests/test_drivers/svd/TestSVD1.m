classdef TestSVD1 < TestSVDecomposer
    
    methods
        function obj = TestSVD1(psd)
            obj = obj@TestSVDecomposer(psd);
        end
        
        function [] = test_fr(obj)
            num_passes = 3;
            block_size = 5;
            logging = 0;
            alg = @(A_, k_, tol_, over_, seed_) svd1(A_, k_, tol_, over_,...
                num_passes, block_size, seed_, logging);

            rank = 15;
            seed = 99;
            alg_tol = 1e-8;
            % Tall matrix
            [U, s, Vt] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            sth = SvdTestHelper.convert(U * diag(S) * Vt, U, s, Vt);

            obj.run_batch(sth, alg, rank - 10, NaN, alg_tol, 0);
            obj.run_batch(sth, alg, rank - 10, NaN, alg_tol, 5);
            obj.run_batch(sth, alg, rank - 3, NaN, alg_tol, 5);

            % Wide matrix
            [U, s, Vt] = gen_test_mat_factors(50, 200, rank, NaN, seed);
            sth = SvdTestHelper.convert(U * diag(S) * Vt, U, s, Vt);

            obj.run_batch(sth, alg, rank - 10, NaN, alg_tol, 0);
            obj.run_batch(sth, alg, rank - 10, NaN, alg_tol, 5);
            obj.run_batch(sth, alg, rank - 3, NaN, alg_tol, 5);
        end

        function [] = test_fp_inexact(obj)
            num_passes = 2;
            block_size = 5;
            alg = @(A_, k_, tol_, over_, seed_) svd1(A_, k_, tol_, over_,...
                num_passes, block_size, seed_, logging);
            seed = 11111;
            rel_err = 0.001;
            rank = 50;

            [U, s, Vt] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            sth = SvdTestHelper.convert(U * diag(S) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2);

            [U, s, ~] = gen_test_mat_factors(50, 200, rank, NaN, seed);
            sth = EigTestHelper.convert(U, s, obj.PSD, seed);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2);
        end

        function [] = test_fp_exact(obj)
            num_passes = 2;
            block_size = 5;
            alg = @(A_, k_, tol_, over_, seed_) svd1(A_, k_, tol_, over_,...
                num_passes, block_size, seed_, logging);
            seed = 1512;
            rank = 15;
            rel_err = 1e-12;

            [U, s, Vt] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            sth = SvdTestHelper.convert(U * diag(S) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2);

            [U, s, Vt] = gen_test_mat_factors(50, 200, rank, NaN, seed);
            sth = SvdTestHelper.convert(U * diag(S) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2);

            [U, s, Vt] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            sth = SvdTestHelper.convert(U * diag(S) * Vt, U, s, Vt);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 0);
            obj.run_batch(sth, alg, rank, rel_err, 1e-8, 2);
        end
    end
end