classdef QbTestHelper
%{
    Objects of this class hold values for algorithm input and output 
    parameters; functions represent unit tests. 
%}
    
    properties
        % Initial data matrix (for tests, randomly generated).
        A 

        % Orthonormal matrix Q of the QB decomposition.
        Q

        % Matrix B of the QB decomposition.
        B

        % Rank of the QB approximation.
        k
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % Initialization function, sets Q and B to empty matrices.
        function obj = QbTestHelper(A, k)
            obj.A = A;
            obj.Q = [];
            obj.B = [];
            obj.k = k;
        end
    
        % Checks if the difference between the initial matrix and
        % approximation is within the specified tolerance.
        function[] = test_exact(obj, tol)
            A = obj.A;
            Q = obj.Q;
            B = obj.B;

            delta = A - Q * B;
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        end
        
        % Checks whether the matrix Q is orthonormal (difference up to a 
        % specified tolerance).
        function[] = test_valid_onb(obj, tol)
            Q = obj.Q;

            gram = Q' * Q;
            delta = gram - eye(size(gram));
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        end
        
        % Checks the quality of computation of matrix B.
        function[] = test_exact_B(obj, tol)
            A = obj.A;
            Q = obj.Q;
            B = obj.B;

            delta = B - Q' * A;
            nrm = norm(delta, 'fro');
            assert(nrm <= tol);
        end
        
        % Check if matrix B has proper rank. 
        function[] = test_exact_rank_B(obj)
            B = obj.B;
            assert(rank(B) == min(size(B))); 
        end
    end
end