classdef TestOSIDecomposer
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function obj = TestOSIDecomposer()
        end

        function[] = run_osid_test(obj, osid, m, n, k, over, p, axis, seed, test_tol, logging)
            
            A = gen_test_mat(m, n, k, 0, 10); 
            [M, P, ~] = osid(A, k, over, p, axis, seed, logging);

            if axis == 0
                A_id = M * A(P, :);
                permuted_coeffs = M(P, :);
                delta_norm = norm(permuted_coeffs - eye(size(P, 1)), 'fro');
                assert(delta_norm < 1e-8)
            elseif axis == 1
                A_id = A(:, P) * M;
                permuted_coeffs = M(:, P);
                delta_norm = norm(permuted_coeffs - eye(size(P, 1)), 'fro');
                assert(delta_norm < 1e-8)
            end
        
            err_rand = norm(A - A_id, 'fro');
            if test_tol < 1e-8
                rel_err = err_rand / norm(A, 'fro');
            else
                [A_id_ref, ~, ~] = obj.reference_osid(A, k, axis);
                err_ref = norm(A - A_id_ref, 'fro');
                rel_err = (err_rand - err_ref) / norm(A, 'fro');
            end
            assert(rel_err < test_tol);
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