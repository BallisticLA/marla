
% Comment: how can i loop through a tructure of objects?
classdef TestQB1 < TestQBecomposer
    
    
    methods
        function obj = TestQB1()
            obj = obj@TestQBecomposer();
        end
        
        function [] = run_var_rank(obj)
            
            spectrum = 0;
            A = [];
            alg_tol = 1e-8;
            test_tol = 1e-8;
            p = 1;
            block_size = 5;
            logging.depth = 0;
            logging.span = 0;

            qth = QbTestHelper(A, 0);
            
            % QB versions.
            alg_1 = @(A, k, seed, logging) rand_qb(A, k, p, seed, logging);
            alg_2 = @(A, k, seed, logging) rand_qb_b(A, block_size, alg_tol, k, p, seed, logging);
            alg_3 = @(A, k, seed, logging) rand_qb_b_pe(A, block_size, alg_tol, k, p, seed, logging);
            alg_4 = @(A, k, seed, logging) rand_qb_sp(A, k, seed, logging);

            % Matrix generator versions.
            gen_1 = @(m, n, k, seed) gen_test_mat(m, n, k, spectrum, seed);
            
            % QB1
            % Cases with low-rank factorization.
            obj.run_exact(qth, alg_1, gen_1, 200, 50, 15, alg_tol, test_tol, logging);
            obj.run_exact(qth, alg_1, gen_1, 50, 200, 15, alg_tol, test_tol, logging);  
            % Cases with full-rank factorizations.
            obj.run_exact(qth, alg_1, gen_1, 200, 50, 50, alg_tol, test_tol, logging); 
            obj.run_exact(qth, alg_1, gen_1, 50, 200, 50, alg_tol, test_tol, logging);     
            % QB2
            % Cases with low-rank factorization.
            obj.run_exact(qth, alg_2, gen_1, 200, 50, 15, alg_tol, test_tol, logging);
            obj.run_exact(qth, alg_2, gen_1, 50, 200, 15, alg_tol, test_tol, logging);  
            % Cases with full-rank factorizations.
            obj.run_exact(qth, alg_2, gen_1, 200, 50, 50, alg_tol, test_tol, logging); 
            obj.run_exact(qth, alg_2, gen_1, 50, 200, 50, alg_tol, test_tol, logging); 
            % QB3
            % Cases with low-rank factorization.
            obj.run_exact(qth, alg_3, gen_1, 200, 50, 15, alg_tol, test_tol, logging);
            obj.run_exact(qth, alg_3, gen_1, 50, 200, 15, alg_tol, test_tol, logging);  
            % Cases with full-rank factorizations.
            obj.run_exact(qth, alg_3, gen_1, 200, 50, 50, alg_tol, test_tol, logging); 
            obj.run_exact(qth, alg_3, gen_1, 50, 200, 50, alg_tol, test_tol, logging); 
            % QB4
            % Cases with low-rank factorization.
            obj.run_exact(qth, alg_4, gen_1, 200, 50, 15, alg_tol, test_tol, logging);
            obj.run_exact(qth, alg_4, gen_1, 50, 200, 15, alg_tol, test_tol, logging);  
            % Cases with full-rank factorizations.
            obj.run_exact(qth, alg_4, gen_1, 200, 50, 50, alg_tol, test_tol, logging); 
            obj.run_exact(qth, alg_4, gen_1, 50, 200, 50, alg_tol, test_tol, logging); 
        end

        function [] = run_overestimated(obj)
            
            spectrum = 0;
            t = 2; % 150
            A = [];
            Q = [];
            B = [];
            alg_tol = 1e-8;
            test_tol = 1e-8;
            p = 1;
            block_size = 5;
            logging.depth = 0;
            logging.span = 0;

            qth = QbTestHelper(A, 0);
            
            % QB versions.
            alg_1 = @(A, k, seed, logging) rand_qb(A, k, p, seed, logging);
            alg_2 = @(A, k, seed, logging) rand_qb_b(A, block_size, alg_tol, k, p, seed, logging);
            alg_3 = @(A, k, seed, logging) rand_qb_b_pe(A, block_size, alg_tol, k, p, seed, logging);
            alg_4 = @(A, k, seed, logging) rand_qb_sp(A, k, seed, logging);

            % Matrix generator versions.
            gen_1 = @(m, n, k, seed) gen_test_mat(m, n, k, spectrum, s);
    
            % Tests with overestimated rank.
            % QB1
            obj.run_estimated(qth, alg_1, gen_1, 200, 50, floor(15 * 1.2), alg_tol, test_tol, logging);
            obj.run_estimated(qth, alg_1, gen_1, 50, 200, floor(15 * 1.2), alg_tol, test_tol, logging);
            % QB2
            obj.run_estimated(qth, alg_2, gen_1, 200, 50, floor(15 * 1.2), alg_tol, test_tol, logging);
            obj.run_estimated(qth, alg_2, gen_1, 50, 200, floor(15 * 1.2), alg_tol, test_tol, logging);
            % QB3
            obj.run_estimated(qth, alg_3, gen_1, 200, 50, floor(15 * 1.2), alg_tol, test_tol, logging);
            obj.run_estimated(qth, alg_3, gen_1, 50, 200, floor(15 * 1.2), alg_tol, test_tol, logging);
            % QB4
            obj.run_estimated(qth, alg_4, gen_1, 200, 50, floor(15 * 1.2), alg_tol, test_tol, logging);
            obj.run_estimated(qth, alg_4, gen_1, 50, 200, floor(15 * 1.2), alg_tol, test_tol, logging);
        end

        function [] = run_fast_decay(obj)
            
            spectrum = 0;
            % Controls the "spped" of singular value decay.
            t = 2;
            A = [];
            Q = [];
            B = [];
            alg_tol = 1e-8;
            test_tol = 1e-8;
            p = 1;
            block_size = 5;
            logging.depth = 0;
            logging.span = 0;

            qth = QbTestHelper(A, 0);
            
            % QB versions.
            alg_1 = @(A, k, seed, logging) rand_qb(A, k, p, seed, logging);
            alg_2 = @(A, k, seed, logging) rand_qb_b(A, block_size, alg_tol, k, p, seed, logging);
            alg_3 = @(A, k, seed, logging) rand_qb_b_pe(A, block_size, alg_tol, k, p, seed, logging);
            alg_4 = @(A, k, seed, logging) rand_qb_sp(A, k, seed, logging);

            % Matrix generator versions.
            gen_2 = @(m, n, k, seed) gen_exp_spectrum(m, n, k, t, seed);

            % QB1
            obj.run_exact(qth, alg_1, gen_2, 200, 50, 15, alg_tol, test_tol, logging); 
            % QB2
            obj.run_exact(qth, alg_2, gen_2, 200, 50, 15, alg_tol, test_tol, logging); 
            % QB3
            obj.run_exact(qth, alg_3, gen_2, 200, 50, 15, alg_tol, test_tol, logging); 
            % QB4
            obj.run_exact(qth, alg_4, gen_2, 200, 50, 15, alg_tol, test_tol, logging); 
        end

        function [] = run_slow_decay(obj)
            
            spectrum = 0;
            % Controls the "spped" of singular value decay.
            t = 150;
            seed = 24;
            A = [];
            Q = [];
            B = [];
            alg_tol = 1e-8;
            test_tol = 1e-8;
            p = 1;
            block_size = 5;
            logging.depth = 0;
            logging.span = 0;

            qth = QbTestHelper(A, 0);
            
            % QB versions.
            alg_1 = @(A, k, seed, logging) rand_qb(A, k, p, seed, logging);
            alg_2 = @(A, k, seed, logging) rand_qb_b(A, block_size, alg_tol, k, p, seed, logging);
            alg_3 = @(A, k, seed, logging) rand_qb_b_pe(A, block_size, alg_tol, k, p, seed, logging);
            alg_4 = @(A, k, seed, logging) rand_qb_sp(A, k, seed, logging);

            % Matrix generator versions.
            gen_2 = @(m, n, k, seed) gen_exp_spectrum(m, n, k, t, seed);
    
            % QB1
            obj.run_exact(qth, alg_1, gen_2, 200, 50, 15, alg_tol, test_tol, logging); 
            % QB2
            obj.run_exact(qth, alg_2, gen_2, 200, 50, 15, alg_tol, test_tol, logging); 
            % QB3
            obj.run_exact(qth, alg_3, gen_2, 200, 50, 15, alg_tol, test_tol, logging); 
            % QB4
            obj.run_exact(qth, alg_4, gen_2, 200, 50, 15, alg_tol, test_tol, logging); 
        end

        function [] = run_s_shaped_decay(obj)
            
            A = [];
            Q = [];
            B = [];
            alg_tol = 1e-8;
            test_tol = 1e-8;
            p = 1;
            block_size = 5;
            logging.depth = 0;
            logging.span = 0;

            qth = QbTestHelper(A, 0);
            
            % QB versions.
            alg_1 = @(A, k, seed, logging) rand_qb(A, k, p, seed, logging);
            alg_2 = @(A, k, seed, logging) rand_qb_b(A, block_size, alg_tol, k, p, seed, logging);
            alg_3 = @(A, k, seed, logging) rand_qb_b_pe(A, block_size, alg_tol, k, p, seed, logging);
            alg_4 = @(A, k, seed, logging) rand_qb_sp(A, k, seed, logging);

            % Matrix generator versions.
            gen_3 = @(m, n, k, seed) gen_s_shaped_spectrum(m, n, k, seed);
    
            % QB1
            obj.run_exact(qth, alg_1, gen_3, 200, 50, 15, alg_tol, test_tol, logging); 
            % QB2
            obj.run_exact(qth, alg_2, gen_3, 200, 50, 15, alg_tol, test_tol, logging); 
            % QB3
            obj.run_exact(qth, alg_3, gen_3, 200, 50, 15, alg_tol, test_tol, logging); 
            % QB4
            obj.run_exact(qth, alg_4, gen_3, 200, 50, 15, alg_tol, test_tol, logging);   
        end
    end
end