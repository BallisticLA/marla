classdef TestTSIDecomposer
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function obj = TestTSIDecomposer()
        end

        function[] = run_tsid_test(obj, tsid, m, n, k, rank, over, p, test_tol, seed, logging)
            
            A = gen_test_mat(m, n, rank, 0, 1); 
            [Z, I, X, J, ~] = tsid(A, k, over, p, seed, logging);

            A_id = Z * (A(I, J) * X);
            err_rand = norm(A - A_id, 'fro');
            if test_tol < 1e-8
                rel_err = err_rand / norm(A, 'fro');
            else
                k = max(size(I));  % should be size(I, 1), but let's be safe.
                [A_id_ref, ~, ~] = obj.reference_osid(A, k, 1);
                err_ref = norm(A - A_id_ref, 'fro');
                rel_err = (err_rand - err_ref) / norm(A, 'fro');
            end
        end

        function[A_id, M, P] = reference_osid(obj, A, k, axis)
            [M, P] = qrcp_osid(A, k, axis);
            if axis == 0
                A_id = M * A(P, :);
            else
                A_id = A(:, P) * M;
            end
        end
    end
end