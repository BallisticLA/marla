function [] = test_cur()
    test_simple_exact()
    test_simple_approx()
end


function [] = test_simple_exact()
    seed = 654;
    p = 1;
    alg = @(A_, k_, s_, seed_) cur1(A_, k_, s_, p, seed_);

    m = 100;
    n = 30;
    run_cur_test(alg, m, n, 24, 25, 3, 1e-12, seed);
    run_cur_test(alg, n, m, 24, 25, 3, 1e-12, seed);

    run_cur_test(alg, m, n, 5, 5, 1, 1e-12, seed);
    run_cur_test(alg, n, m, 5, 5, 1, 1e-12, seed);

end


function [] = test_simple_approx()
    seed = 192;
    p = 1;
    alg = @(A_, k_, s_, seed_) cur1(A_, k_, s_, p, seed_);

    m = 100;
    n = 30;
    run_cur_test(alg, m, n, 30, 27, 3, 1e-1, seed);
    run_cur_test(alg, n, m, 30, 27, 3, 1e-1, seed);

    run_cur_test(alg, m, n, 30, 25, 4, 1e-1, seed);
    run_cur_test(alg, n, m, 30, 25, 4, 1e-1, seed);

    run_cur_test(alg, m, n, 30, 5, 5, 3e-1, seed);
    run_cur_test(alg, n, m, 30, 5, 5, 3e-1, seed);    
end


function [] = run_cur_test(alg, m, n, rank, k, over, test_tol, seed)
    A = gen_test_mat(m, n, rank, NaN, seed);
    [Js, U, Is] = alg(A, k, over, seed);
    A_id = A(:, Js) * U * A(Is, :);

    err_rand = norm(A - A_id, 'fro');
    if test_tol < 1e-8
        rel_err = err_rand / norm(A, 'fro');
    else
        [A_id_ref, ~, ~] = reference_osid(A, k, 0);
        err_ref = norm(A - A_id_ref, 'fro');
        rel_err = (err_rand - err_ref) / norm(A, 'fro');
    end
    disp(rel_err)
    assert(rel_err < test_tol)
end

function[A_id, M, P] = reference_osid(A, k, axis)
    addpath('../../comps/interpolative/');
    [M, P] = qrcp_osid(A, k, axis);
    if axis == 0
        A_id = M * A(P, :);
    else
        A_id = A(:, P) * M;
    end
end