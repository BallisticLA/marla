classdef TestTSIDs < TestTSIDecomposer

    methods
        function obj = TestTSIDs()
            obj = obj@TestTSIDecomposer();
        end

        function [] = test_simple_exact(obj)
            alg = @(A, k, over, p, s, logging) tsid1(A, k, over, p, s, logging);
            
            seed = 1;
            test_tol = 1e-12;
            logging.depth = 0;
            logging.span = 0;
    
            % Proceed with a loop when more algorithms added.
            obj.run_tsid_test(alg, 1000, 300, 290, 300, 0, 1, test_tol, seed, logging);
            % ^ Once we add more tests (e.g., that A[Is,Js]
            % is nonsingular), this case might fail.
            obj.run_tsid_test(alg, 100, 30, 30, 30, 0, 1, test_tol, seed, logging);
        end
    
        function [] = test_simple_approx(obj)
            alg = @(A, k, over, p, s, logging) tsid1(A, k, over, p, s, logging);
            
            seed = 2;
            test_tol = 0.05;
            logging.depth = 0;
            logging.span = 0;
            
            % Proceed with a loop when more algorithms added.
            % Tall matrices.
            obj.run_tsid_test(alg, 100, 30, 27, 30, 3, 2, test_tol, seed, logging);  
            obj.run_tsid_test(alg, 100, 30, 25, 30, 4, 2, test_tol, seed, logging);
            
            % Wide matrices.
            obj.run_tsid_test(alg, 30, 100, 27, 30, 3, 2, test_tol, seed, logging);
            obj.run_tsid_test(alg, 30, 100, 25, 30, 4, 2, test_tol, seed, logging);
        end
    end
end