classdef TestCURs < TestCURecomposer

    methods
        function obj = TestCURs()
            obj = obj@TestCURecomposer();
        end

        function [] = test_simple_exact(obj)
            alg = @(A, k, over, p, s, logging) cur1(A, k, over, p, s, logging);
            
            seed = 654;
            test_tol = 1e-12;
            % Test matrix dimensions.
            m = 100;
            n = 30;
            logging.depth = 0;
            logging.span = 0;
    
            obj.run_cur_test(alg, m, n, 24, 25, 3, 1, test_tol, seed, logging);
            obj.run_cur_test(alg, n, m, 24, 25, 3, 1, test_tol, seed, logging);
        
            obj.run_cur_test(alg, m, n, 5, 5, 1, 1, test_tol, seed, logging);
            obj.run_cur_test(alg, n, m, 5, 5, 1, 1, test_tol, seed, logging);
        end
    
        function [] = test_simple_approx(obj)
            alg = @(A, k, over, p, s, logging) cur1(A, k, over, p, s, logging);
            
            seed = 192;
            % Test matrix dimensions.
            m = 100;
            n = 30;
            logging.depth = 0;
            logging.span = 0;
    
            obj.run_cur_test(alg, m, n, 30, 27, 3, 1, 1e-1, seed, logging);
            obj.run_cur_test(alg, n, m, 30, 27, 3, 1, 1e-1, seed, logging);
        
            obj.run_cur_test(alg, m, n, 30, 25, 4, 1, 1e-1, seed, logging);
            obj.run_cur_test(alg, n, m, 30, 25, 4, 1, 1e-1, seed, logging);
        
            obj.run_cur_test(alg, m, n, 30, 5, 5, 1, 3e-1, seed, logging);
            obj.run_cur_test(alg, n, m, 30, 5, 5, 1, 3e-1, seed, logging);         
        end
    end
end