classdef TestOSIDs < TestOSIDecomposer

    methods
        function obj = TestOSIDs()
            obj = obj@TestOSIDecomposer();
        end

        function [] = test_simple_exact(obj)
            alg_1 = @(A, k, over, p, axis, seed, logging) osid1(A, k, over, p, axis, seed, logging);
            alg_2 = @(A, k, over, p, axis, seed, logging) osid2(A, k, over, p, axis, seed, logging);
            
            seed = 487346485;
            test_tol = 1e-12;
            % Test matrix dimensions.
            m = 100;
            n = 30;
            logging.depth = 0;
            logging.span = 0;
    
            for i = 1 : 2
                if i == 1, alg = alg_1; else, alg = alg_2; end
                % Row IDs
                obj.run_osid_test(alg, m, n, 29, 0, 1, 0, seed, test_tol, logging);
                obj.run_osid_test(alg, m, n, 28, 1, 1, 0, seed, test_tol, logging);
                obj.run_osid_test(alg, m, n, 3, 0, 1, 0, seed, test_tol, logging);
                % Column IDs
                obj.run_osid_test(alg, m, n, 29, 0, 1, 1, seed, test_tol, logging);
                obj.run_osid_test(alg, m, n, 28, 1, 1, 1, seed, test_tol, logging);
                obj.run_osid_test(alg, m, n, 3, 0, 1, 1, seed, test_tol, logging);
            end
        end
    
        function [] = test_simple_approx(obj)
            alg_1 = @(A, k, over, p, axis, seed, logging) osid1(A, k, over, p, axis, seed, logging);
            alg_2 = @(A, k, over, p, axis, seed, logging) osid2(A, k, over, p, axis, seed, logging);
            
            seed = 93;
            test_tol = 0.1;
            % Test matrix dimensions.
            m = 100;
            n = 30;
            logging.depth = 0;
            logging.span = 0;
    
            for i = 1 : 2
                if i == 1, alg = alg_1; else, alg = alg_2; end
                % Tall matrices.
                obj.run_osid_test(alg, m, n, 27, 3, 1, 1, seed, test_tol, logging);
                obj.run_osid_test(alg, m, n, 25, 4, 1, 1, seed, test_tol, logging);
                obj.run_osid_test(alg, m, n, 5, 5, 1, 0, seed, test_tol, logging);
                
                % Wide matrices.
                obj.run_osid_test(alg, n, m, 27, 3, 1, 1, seed, test_tol, logging);
                obj.run_osid_test(alg, n, m, 25, 4, 1, 1, seed, test_tol, logging);
                obj.run_osid_test(alg, n, m, 5, 5, 1, 0, seed, test_tol, logging);
            end        
        end
    end
end