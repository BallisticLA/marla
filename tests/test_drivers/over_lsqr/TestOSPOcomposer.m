classdef TestOSPOcomposer
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function obj = TestOSPOcomposer()
        end
        function[] = run_inconsistent(obj, ath, spo, alg_tol,...
                         test_tol)

            % Call the algorithms
            [x_approx, ~] = spo(ath.A, ath.b);
            ath.x_approx = x_approx;

            % Run tests
            ath.test_residual_proj(test_tol);
            ath.test_x_angle(test_tol);
            ath.test_x_norm(test_tol);
        end

        function[] = run_consistent(obj, ath, spo, alg_tol,...
                         test_tol)

            % Call the algorithms
            [x_approx, ~] = spo(ath.A, ath.b);
            ath.x_approx = x_approx;

            % Run tests
            ath.test_x_angle(test_tol);
            ath.test_x_norm(test_tol);
            ath.test_objective(test_tol);
        end
    end
end