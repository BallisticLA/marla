classdef SpoTestHelper
    
    properties
        A
        b
        x_opt
        x_approx
        U
        s
        Vt
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function self = SpoTestHelper(A, b, x_opt, U, s, Vt)
            self.A = A;
            self.b = b;
            self.x_opt = x_opt;
            self.x_approx = [];
            self.U = U;
            self.s = s;
            self.Vt = Vt;
        end
        
        function[] = test_x_angle(self, tol)
            % x' x_opt >= (1 - tol)*||x|| ||x_opt||
            y_opt = self.Vt * self.x_opt;
            norm_y_opt = norm(y_opt, 2);
            y = self.Vt * self.x_approx;
            norm_y = norm(y, 2);
            if norm_y_opt < 1e-8
                % Norm is too small to accurately compute cosine
                assert(abs(norm_y - norm_y_opt) <= tol);
            else
                y_opt = y_opt / norm_y_opt;
                y = y / norm_y;
                cosine = dot(y, y_opt);
                assert(cosine >= (1 - tol));
            end
        end

        function[] = test_x_norm(self, tol)
            % (1 - tol)*||x_opt|| <= ||x|| <= (1+tol)*||x_opt|| + tol
            nrm = norm(self.Vt * self.x_approx, 2);
            norm_opt = norm(self.Vt * self.x_approx, 2);
            assert(nrm <= ((1+tol)*norm_opt + tol));
            assert(((1-tol)*norm_opt) <= nrm);    
        end
        
        function[] = test_delta_x(self, tol)
            %||x - x_opt|| <= tol
            delta_x = self.x_opt - self.x_approx;
            nrm = norm(delta_x, 2) / (1 + min(norm(self.x_opt, 2), norm(self.x_approx, 2)));
            assert(nrm <= tol);
        end
        
        function[] = test_residual_proj(self, tol)
            % || U U' (A x - b) || / ||A x - b|| <= tol
            % This test is probably better scaled than the normal equations
            residual = self.A * self.x_approx - self.b;
            residual_proj = self.U * (self.U' * residual);
            nrm = norm(residual_proj, 2) / norm(residual, 2);
            assert(nrm <= tol);
        end
        
        function[] = test_objective(self, tol)
            % ||A x - b|| <= ||A x_opt - b|| + tol
            res_approx = self.b - self.A * self.x_approx;
            res_opt = self.b - self.A * self.x_opt;
            nrm_approx = norm(res_approx, 2);
            nrm_opt = norm(res_opt, 2);
            assert(nrm_approx <= (nrm_opt + tol));
        end
        
        function[] = test_normal_eqs(self, tol)
            % || A' A x - A' b|| <= tol
            gap = (self.A)' * self.b - (self.A)' * (self.A * self.x_approx);
            nrm = norm(gap, 2);
            assert(nrm <= tol);
        end
    end
end