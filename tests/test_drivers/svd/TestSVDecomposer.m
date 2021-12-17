classdef TestSVDecomposer
    
    properties
        SEEDS = [38972, 653, 1222];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function obj = TestSVDecomposer()
        end

        function [] = run_batch(obj, sth, rsvd, target_rank,...
                target_tol, test_tol, over)
            % obj: TestSVDecomposer
            % sth: SvdTestHelper
            % rsvd: (partial) function handle for svd1
            for idx = 1:length(obj.SEEDS)
                seed = obj.SEEDS(idx);

                % Call the algorithm
                [U, s, Vt, ~] = rsvd(sth.A,...
                    target_rank, target_tol, over, seed);
                sth.UsVt = [U, s, Vt];
                
                % Run tests
                sth.test_conformable();
                sth.test_valid_singvals();
                sth.test_valid_onb(test_tol);
                test_abs_fro_error(target_tol);
            end
        end
    end

end