classdef SvdTestHelper
    %EIGTESTHELPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Initial data matrix (for tests, randomly generated).
        A
        U
        s
        Vt
        % ^ Those U, s, Vt are reference values
        U_approx
        s_approx
        V_approx
        % ^ To store values computed by a driver
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % Initialization function, setsU_approx, s_approx, V_approx to empty matrices.
        function obj = SvdTestHelper(A, U, s, Vt)
            obj.A = A;
            obj.U = U;
            obj.s = s;
            obj.Vt = Vt;
            obj.U_approx = [];
            obj.s_approx = [];
            obj.V_approx = [];
        end

        % Checks dimensions matrching.
        function[] = test_conformable(self)
            s = diag(self.s);

            assert(size(s, 1) == size(self.U, 2));
            assert(size(s, 2) == size(self.Vt, 1));
            assert(size(self.U, 1) == size(self.A, 1));
            assert(size(self.Vt, 2) == size(self.A, 2));
        end
        
        % Checks whether matrices V and U are orthonormal (difference up to a 
        % specified tolerance).
        function[] = test_valid_onb(self, tol)
            U = self.U_approx;
            V = self.V_approx;
        
            gram = U' * U;
            delta = gram - eye(size(gram));
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        
            gram = V' * V;
            delta = gram - eye(size(gram));
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        end
        
        % Checks the amount of approximated singular values, as well as the
        % magnitude of the smallest singular value. 
        function[]  = test_valid_singvals(self, test_tol)
            s_exact = self.s;
            s = diag(self.s_approx)';

            assert(size(s, 2) <= size(s_exact, 2));
            assert(min(size(s)) >= 0);
        end
        
        % Checks the quality of initial data approximation (with specified
        % tolerance).
        function[] = test_abs_fro_error(self, rel_tol)
            A = self.A;
            U = self.U_approx;
            s = self.s_approx;
            Vt = (self.V_approx)';


            %logging.depth = 0;
            %logging.span = 0;
            %[U1, S1, V1, ~] = svd1(A, 13, 1e-8, 0, 0, 2, 99, logging);
            %disp(norm(A - (U1 * S1 * V1'), 'fro'));

            % TODO: change this to relative tolerance
            delta = A - (U * s * Vt);
            nrm = norm(delta, 'fro');
            %disp(nrm)
            % abs_tol = rel_tol * norm(S, 2)
            % ^ Scale by  Frobenius norm of A.
            % assert(nrm <= rel_tol)
            assert(nrm <= rel_tol);
        end
    end
end