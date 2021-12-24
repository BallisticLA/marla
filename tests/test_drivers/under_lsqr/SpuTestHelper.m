classdef SpuTestHelper
    
    properties
        A
        c
        y_opt
        y_approx
        s
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function self = SpuTestHelper(A, c, y_opt, s)
            self.A = A;
            self.c = c;
            self.y_opt = y_opt;
            self.y_approx = [];
            self.s = s;
        end

        function[] = test_delta_y(self, tol)
            % ||y - y_opt|| <= tol
            delta_y = self.y_opt - self.y_approx;
            nrm = norm(delta_y, 2) / (1 + min(norm(self.y_opt, 2), norm(self.y_approx, 2)));
            disp(nrm); 
            assert(nrm <= tol);
        end
        
        function[] = test_objective(self, tol)
            % ||y|| <= ||y_opt|| + tol
            nrm_opt = norm(self.y_opt, 2);
            nrm_approx = norm(self.y_approx, 2);
            assert(nrm_approx <= (1 + tol)*nrm_opt);
        end
        
        function[] = test_residual(self, tol)
            % ||A' y - c|| <= tol
            res = self.A' * self.y_approx - self.c;
            nrm = norm(res, 2);
            assert(nrm <= tol);
        end
    end
end