classdef TestSPU1 < TestSPUcomposer
    
    methods
        function obj = TestSPU1()
            obj = obj@TestSPUcomposer();
        end
        
        function [] = linspace_spec(obj)
            % Test data generation.
            m = 1000;
            n = 100;
            cond_num = 1e5;
            spectrum = linspace(cond_num^0.5, cond_num^(-0.5), n);
            seed = MarlaRandStream(3893874);
            A = gen_test_mat(m, n, n, spectrum, seed);
            y0 = randn(seed, m, 1);
            c = A' * y0;
            ath = SpuTestHelper(A, c, A' \ c, spectrum);

            % Parameters.
            seed = 998765;
            sampling_factor = 3;
            alg_tol = 1e-12;
            iter_lim = 50;
            logging.depth = 0;
            logging.span = 0;

            alg = @(A, b) spu1(A, c, sampling_factor, alg_tol, iter_lim, logging, seed);
            obj.run_test(ath, alg, alg_tol, 1e-6);
        end

        function [] = logspace_spec(obj)
            % Test data generation.
            m = 1000;
            n = 100;
            cond_num = 1e5;
            spectrum = logspace(log10(cond_num)/2, -log10(cond_num)/2, n);
            seed = MarlaRandStream(3893874);
            A = gen_test_mat(m, n, n, spectrum, seed);
            y0 = randn(seed, m, 1);
            c = A' * y0;
            ath = SpuTestHelper(A, c, A' \ c, spectrum);

            % Parameters.
            seed = 998765;
            sampling_factor = 3;
            alg_tol = 1e-12;
            iter_lim = 50;
            logging.depth = 0;
            logging.span = 0;

            alg = @(A, b) spu1(A, c, sampling_factor, alg_tol, iter_lim, logging, seed);
            obj.run_test(ath, alg, alg_tol, 1e-6);
        end

        function [] = lowrank_linspace_spec(obj)
            % Test data generation.
            m = 1000;
            n = 100;
            cond_num = 1e5;
            rk = 80;
            spectrum = linspace(cond_num^0.5, cond_num^(-0.5), rk);
            seed = MarlaRandStream(3893874);
            A = gen_test_mat(m, n, rk, spectrum, seed);
            y0 = randn(seed, m,1);
            c = A' * y0;
            ath = SpuTestHelper(A, c, A' \ c, spectrum);

            % Parameters.
            seed = 998765;
            sampling_factor = 3;
            alg_tol = 1e-12;
            iter_lim = 50;
            logging.depth = 0;
            logging.span = 0;

            alg = @(A, b) spu1(A, c, sampling_factor, alg_tol, iter_lim, logging, seed);
            obj.run_test(ath, alg, alg_tol, 1e-6);
        end
    end
end