classdef EigTestHelper
    %EIGTESTHELPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A
        V
        lamb
        fro_A
        V_approx
        lamb_approx
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function obj = EigTestHelper(A, V, lambda)
            obj.A = A;
            obj.fro_A = norm(A, 'fro');
            obj.V = V;
            obj.lamb = lambda;
        end
        
        function [] = test_conformable(obj)
            assert(size(obj.lamb_approx,1) == size(obj.V_approx, 2));
            assert(size(obj.V_approx, 1) == size(obj.A, 1));
        end

        function [] = test_valid_onb(obj, fro_tol)
            V = obj.V_approx;
            lamb = obj.lamb_approx;
            gram_V = V' * V;
            delta_V = gram_V - eye(size(lamb, 1));
            nrm_V = norm(delta_V, 'fro');
            assert(nrm_V <= fro_tol);
        end

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

