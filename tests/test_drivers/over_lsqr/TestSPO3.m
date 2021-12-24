classdef TestSPO3 < TestOSPOcomposer
    
    methods
        function obj = TestSPO3()
            obj = obj@TestOSPOcomposer();
        end
        
        function [] = consistent_tall(obj)
            % Test data generation.
            stream = MarlaRandStream(98743);
            A = gen_simp_mat(100, 10, 1, stream);
            [U, s, Vt] = svd(A);
            x = randn(10, 1);      
            ath = SpoTestHelper(A, A * x, x, U, s, Vt);

            % Parameters.
            seed = 0;
            sampling_factor = 1;
            alg_tol = 0.0;
            iter_lim = 1;
            use_chol = false;
            logging.depth = 0;
            logging.span = 0;

            alg = @(A, b) spo3(A, b, sampling_factor, alg_tol, iter_lim, use_chol, logging, seed);
            obj.run_consistent(ath, alg, alg_tol, 1e-12);
    
            use_chol = true;

            alg = @(A, b) spo3(A, b, sampling_factor, alg_tol, iter_lim, use_chol, logging, seed);
            obj.run_consistent(ath, alg, alg_tol, 1e-12);
        end

        function [] = consistent_square(obj)
            % Test data generation.
            stream = MarlaRandStream(98743);
            m = 100;
            n = 10;
            s = diag(rand(stream, n, 1) + 1e-4);
            U = randn(stream, m, n);
            [U,~] = qr(U, 0);
            V = randn(stream, n, n);
            [V, ~] = qr(V, 0);
            Vt = V'; 
            A = (U * s) * Vt;
            x = randn(stream, n, 1);
            ath = SpoTestHelper(A, A * x, x, U, s, Vt);

            % Parameters.
            seed = 0;
            sampling_factor = 1;
            alg_tol = 0.0;
            iter_lim = 1;
            use_chol = false;
            logging.depth = 0;
            logging.span = 0;

            alg = @(A, b) spo3(A, b, sampling_factor, alg_tol, iter_lim, use_chol, logging, seed);
            obj.run_consistent(ath, alg, alg_tol, 1e-12);

            use_chol = true;

            alg = @(A, b) spo3(A, b, sampling_factor, alg_tol, iter_lim, use_chol, logging, seed);
            obj.run_consistent(ath, alg, alg_tol, 1e-12);
        end

        function [] = inconsistent_orthog(obj)
            % Test data generation.
            stream = MarlaRandStream(98743);
            m = 1000;
            n = 100;
            U = randn(stream, m, n);
            [U,~] = qr(U, 0);
            V = randn(stream, n, n);
            [V, ~] = qr(V, 0);
            Vt = V'; 
            s = diag(rand(stream, n, 1) + 1e-4);
            A = U * s * Vt;
            b = randn(stream, m, 1);
            b = b - U * (U' * b);
            ath = SpoTestHelper(A, b * 1e2 / norm(b, 2), zeros(n, 1), U, s, Vt);

            % Parameters.
            seed = 0;
            sampling_factor = 3;
            alg_tol = 1e-12;
            iter_lim = 100;
            use_chol = false;
            logging.depth = 0;
            logging.span = 0;

            alg = @(A, b) spo3(A, b, sampling_factor, alg_tol, iter_lim, use_chol, logging, seed);
            obj.run_inconsistent(ath, alg, alg_tol, 1e-6);

            use_chol = true;

            alg = @(A, b) spo3(A, b, sampling_factor, alg_tol, iter_lim, use_chol, logging, seed);
            obj.run_inconsistent(ath, alg, alg_tol, 1e-6);
        end

        function [] = inconsistent_gen(obj)
            % Test data generation.
            stream = MarlaRandStream(98743);
            m = 1000;
            n = 100;
            num_hi = 30;
            num_lo = n - num_hi;
            % Make A
            hi_spec = 1e5 * ones(1, num_hi) + rand(1, num_hi);
            lo_spec = ones(1, num_lo) + rand(1, num_lo);
            spec = diag(cat(2, hi_spec, lo_spec));
            U = randn(stream, m, n);
            [U,~] = qr(U, 0);
            V = randn(stream, n, n);
            [V,~] = qr(V, 0);
            Vt = V';
            A = (U * spec) * Vt;
            % Make b
            hi_x = randn(stream, num_hi, 1) / 1e5;
            lo_x = randn(stream, num_lo, 1);
            x = cat(1, hi_x, lo_x);
            b_orth = randn(stream, m, 1) * 1e2;
            % orthogonal to range(A)
            b_orth = b_orth - U * (U' * b_orth);
            ath = SpoTestHelper(A, A * x + b_orth, x, U, spec, Vt);

            % Parameters.
            seed = 0;
            sampling_factor = 3;
            alg_tol = 1e-12;
            iter_lim = 100;
            use_chol = false;
            logging.depth = 0;
            logging.span = 0;

            alg = @(A, b) spo3(A, b, sampling_factor, alg_tol, iter_lim, use_chol, logging, seed);
            obj.run_inconsistent(ath, alg, alg_tol, 1e-6);

            use_chol = true;

            alg = @(A, b) spo3(A, b, sampling_factor, alg_tol, iter_lim, use_chol, logging, seed);
            obj.run_inconsistent(ath, alg, alg_tol, 1e-6);
        end
    end
end