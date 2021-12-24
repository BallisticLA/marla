classdef TestSPUcomposer
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function obj = TestSPUcomposer()
        end
        function[] = run_test(obj, ath, spo, alg_tol,...
                         test_tol)

            % Call the algorithms
            [y_approx, ~] = spo(ath.A, ath.c);
            ath.y_approx = y_approx;

            % Run tests
            ath.test_objective(test_tol);
            ath.test_residual(test_tol);
            % ath.test_delta_y(test_tol);  % skip this (here, and in Python)
        end
    end
end