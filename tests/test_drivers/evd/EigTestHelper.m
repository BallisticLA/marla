classdef EigTestHelper
%{
    Objects of thsi class hold values for algorithm input and output 
    parameters; functions represent unit tests. 
%}
    
    properties
        % Initial data matrix (for tests, randomly generated).
        A

        % Orthonormal matrix
        V

        % Vector of eigenvalues
        lamb

        % Frobenius norm of A
        fro_A

        % Approximations
        V_approx
        lamb_approx
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        % Initialization function.
        function obj = EigTestHelper(A, V, lambda)
            obj.A = A;
            obj.fro_A = norm(A, 'fro');
            obj.V = V;
            obj.lamb = lambda;
        end
        
        % Checks matching dimensions. 
        function [] = test_conformable(obj)
            assert(size(obj.lamb_approx,1) == size(obj.V_approx, 2));
            assert(size(obj.V_approx, 1) == size(obj.A, 1));
        end

        % Checks whether the matrix V is orthonormal (difference up to a 
        % specified tolerance).
        function [] = test_valid_onb(obj, fro_tol)
            V = obj.V_approx;
            lamb = obj.lamb_approx;
            gram_V = V' * V;
            delta_V = gram_V - eye(size(lamb, 1));
            nrm_V = norm(delta_V, 'fro');
            assert(nrm_V <= fro_tol);
        end
        
        % Checks the size of eigenvalue vector, as well as magnitude of the
        % smallest eigenvalue.
        function [] = test_eigvals(obj, psd)
            lamb = obj.lamb_approx;
            assert(length(lamb) <= length(obj.lamb));
            abslamb_rev = flipud(abs(lamb));
            diffs = diff(abslamb_rev);
            assert(min(diffs) >= 0.0);
            if psd
                assert(min(lamb) >= 0);
            end
        end

        % Checks the approximation quality (up to a specified tolerance).
        function [] = test_fro_error(obj, rel_fro_tol)
            V = obj.V_approx;
            lamb = obj.lamb_approx;
            delta = obj.A - V * diag(lamb) * V';
            nrm = norm(delta, 'fro');
            rel_limit = rel_fro_tol * obj.fro_A;
            assert(nrm <= rel_limit);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Static methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initial data conevsion from SVD to EVD form. 
    methods(Static)
        function eth = convert(U, spectrum, psd, s)
            if psd
                lamb = spectrum;
            else
                s = MarlaRandStream(s);
                flips = rand(s, size(spectrum)) < 0.5;
                signs = ones(size(spectrum));
                signs(flips) = -1;
                lamb = spectrum .* signs;
            end
            A = U * diag(lamb) * U';
            eth = EigTestHelper(A, U, lamb);
        end
    end
end

