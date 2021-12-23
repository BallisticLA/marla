classdef QbTestHelper
    %EIGTESTHELPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A
        Q
        B
        k
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function obj = QbTestHelper(A, k)
            obj.A = A;
            obj.Q = [];
            obj.B = [];
            obj.k = k;
        end
        
        function [] = test_conformable(obj)
            assert(size(obj.lamb_approx,1) == size(obj.V_approx, 2));
            assert(size(obj.V_approx, 1) == size(obj.A, 1));
        end

        function[] = test_exact(obj, tol)
            A = obj.A;
            Q = obj.Q;
            B = obj.B;

            delta = A - Q * B;
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        end
        
        function[] = test_valid_onb(obj, tol)
            Q = obj.Q;

            gram = Q' * Q;
            delta = gram - eye(size(gram));
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        end
        
        function[] = test_exact_B(obj, tol)
            A = obj.A;
            Q = obj.Q;
            B = obj.B;

            delta = B - Q' * A;
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        end
        
        function[] = test_exact_rank_B(obj)
            B = obj.B;
            assert(rank(B) == min(size(B))); 
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Static methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Static)
        function qth = convert(A, Q, B)
            qth = QbTestHelper(A, Q, B);
        end
    end
end