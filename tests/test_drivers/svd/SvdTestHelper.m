classdef SvdTestHelper
    %EIGTESTHELPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A
        U
        s
        Vt
        % ^ Those U, s, Vt are reference values
        UsVt
        % ^ To store values computed by a driver
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function obj = SvdTestHelper(A, U, s, Vt)
            obj.A = A;
            obj.U = U;
            obj.s = s;
            obj.Vt = Vt;
            obj.UsVt = [];
        end

        function[] = test_conformable(self)
            assert(size(self.s, 1) == size(self.U, 2));
            assert(size(self.s, 2) == size(self.Vt, 2));
            assert(size(self.U, 1) == size(self.A, 1));
            assert(size(self.Vt, 1) == size(self.A, 2));
        end
        
        function[] = test_valid_onb(self, tol)
            [U, ~, Vt] = self.UsVt;
        
            gram = U' * U;
            delta = gram - eye(size(gram));
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        
            gram = Vt * Vt';
            delta = gram - eye(size(gram));
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        end
        
        function[]  = test_valid_singvals(self, test_tol)
            self.s_exact = self.s;
            [~, s, ~] = self.UsVt;
            assert(size(self.s, 1) <= size(self.s_exact, 1));
            assert(min(size(self.s)) >= 0);
        end
        
        function[] = test_abs_fro_error(self, rel_tol)
            A = self.A;
            [U, s, Vt] = self.UsVt;
            % TODO: change this to relative tolerance
            delta = A - (U * s) * Vt;
            nrm = norm(delta, 'fro');
            % abs_tol = rel_tol * norm(S, 2)
            % ^ Scale by  Frobenius norm of A.
            % assert(nrm <= rel_tol)
            assert(nrm <= rel_tol);
        end
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Static methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Static)
        function sth = convert(A, U, s, Vt)
            sth = SvdTestHelper(A, U, s, Vt);
        end
    end
end