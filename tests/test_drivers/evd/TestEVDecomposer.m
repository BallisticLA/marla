classdef TestEVDecomposer
    
    properties
        SEEDS = [1,2,3];
        PSD = false;
        INFLATE_TEST_TOL = 1.0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function obj = TestEVDecomposer(psd)
            obj.PSD = psd;
        end

        function [] = run_batch(obj, eth, revd, target_rank,...
                target_tol, test_tol, over, logging)
            % obj: TestEVDecomposer
            % eth: EigTestHelper
            % revd: (partial) function handle for evd1 or evd2
            for idx = 1:length(obj.SEEDS)
                seed = obj.SEEDS(idx);

                % Call the algorithm
                [V_approx, lamb_approx, ~] = revd(eth.A,...
                    target_rank, target_tol, over, seed, logging);
                eth.V_approx = V_approx;
                eth.lamb_approx = lamb_approx;
                
                % Run tests
                eth.test_conformable();
                eth.test_eigvals(obj.PSD);
                eth.test_valid_onb(test_tol);
                if ~isnan(target_tol)
                    reltol = obj.INFLATE_TEST_TOL * target_tol;
                    eth.test_fro_error(reltol);
                end
            end
        end
    end

end

