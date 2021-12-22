classdef TestCURecomposer
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Instance methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function obj = TestCURecomposer()
        end

        function[] = run_cur_test(obj, cur, m, n, rank, k, over, p, test_tol, seed, logging)
            
            A = gen_test_mat(m, n, rank, NaN, seed);
            [Js, U, Is, ~] = cur(A, k, over, p, seed, logging);
            A_id = A(:, Js) * U * A(Is, :);
        
            err_rand = norm(A - A_id, 'fro');
            if test_tol < 1e-8
                rel_err = err_rand / norm(A, 'fro');
            else
                [A_id_ref, ~, ~] = obj.reference_osid(A, k, 0);
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