classdef TestEVD2 < TestEVDecomposer
    
    methods
        function obj = TestEVD2()
            obj = obj@TestEVDecomposer(true);
        end
        
        function [] = test_fr(obj)
            num_passes = 3;
            logging.depth = 0;
            logging.span = 0;
            alg = @(A_, k_, tol_, over_, seed_, logging) evd2(A_, k_, over_,...
                num_passes, seed_, logging);

            seed = 93641;
            rank = 15;
            [U, s, ~] = gen_test_mat_factors(200, 50, rank, NaN, seed);
            eth = EigTestHelper.convert(U, s, obj.PSD, seed);

            obj.run_batch(eth, alg, rank - 10, NaN, 1e-8, 0, logging);
            obj.run_batch(eth, alg, rank - 10, NaN, 1e-8, 5, logging);
            obj.run_batch(eth, alg, rank - 3, NaN, 1e-8, 5, logging);
        end
    end
end
