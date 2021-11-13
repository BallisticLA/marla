function[] = test_osid()
    seed = 10;
    addpath('../matrix_generators/');

    test1.k = 29;
    test1.A = gen_test_mat(100, 30, test1.k, 0, seed); 

    test2.k = 28;
    test2.A = gen_test_mat(100, 30, test2.k, 0, seed); 

    test3.k = 3;
    test3.A = gen_test_mat(100, 30, test3.k, 0, seed); 

    test4.k = 30;
    test4.A = gen_test_mat(100, 30, test4.k, 0, seed); 

    test5.k = 30;
    test5.A = gen_test_mat(30, 100, test5.k, 0, seed); 

    test_simple_exact(test1, test2, test3);
    test_simple_approx(test4, test5);
end

function[] = test_simple_exact(test1, test2, test3)
    seed = 487346485;
    addpath('../../drivers/interpolative/')
    % osid1
    % Row IDs
    [M, P] = osid1(test1.A, 29, 0, 1, 0, seed);
    run_osid_test(test1.A, M, P, 29, 0, 1e-12);

    [M, P] = osid1(test2.A, 28, 1, 1, 0, seed);
    run_osid_test(test2.A, M, P, 28, 0, 1e-12);
    % ^ The test above originally failed because A went into one-sided ID
    % with exact rank 28, then one-sided ID ran a rangefinder
    % which bumped the rank up to 29; since we performed one-sided
    % ID on the orthogonal matrix Q there was no indication that the
    % sketch which led to Q was actually near-singular.
    % Conclusion: have one-sided ID use a row sketcher (as in PARLA)
    % not a rangefinder (which is often used in research papers).

    [M, P] = osid1(test3.A, 3, 0, 1, 0, seed);
    run_osid_test(test3.A, M, P, 3, 0, 1e-12);
    % Column IDs
    [M, P] = osid1(test1.A, 29, 0, 1, 1, seed);
    run_osid_test(test1.A, M, P, 29, 1, 1e-12);

    [M, P] = osid1(test2.A, 28, 1, 1, 1, seed);
    run_osid_test(test2.A, M, P, 28, 1, 1e-12);

    [M, P] = osid1(test3.A, 3, 0, 1, 1, seed);
    run_osid_test(test3.A, M, P, 3, 1, 1e-12);

    % osid2
    % Row IDs
    [M, P] = osid2(test1.A, 29, 0, 1, 0, seed);
    run_osid_test(test1.A, M, P, 29, 0, 1e-12);

    [M, P] = osid2(test2.A, 28, 1, 1, 0, seed);
    run_osid_test(test2.A, M, P, 28, 0, 1e-12);

    [M, P] = osid2(test3.A, 3, 0, 1, 0, seed);
    run_osid_test(test3.A, M, P, 3, 0, 1e-12);
    % Column IDs
    [M, P] = osid2(test1.A, 29, 0, 1, 1, seed);
    run_osid_test(test1.A, M, P, 29, 1, 1e-12);

    [M, P] = osid2(test2.A, 28, 1, 1, 1, seed);
    run_osid_test(test2.A, M, P, 28, 1, 1e-12);

    [M, P] = osid2(test3.A, 3, 0, 1, 1, seed);
    run_osid_test(test3.A, M, P, 3, 1, 1e-12);
end

function[] = test_simple_approx(test4, test5)
    seed = 93;
    addpath('../../drivers/interpolative/')
    % osid1
    % Tall matrices.
    [M, P] = osid1(test4.A, 27, 3, 1, 1, seed);
    run_osid_test(test4.A, M, P, 29, 1, 0.1);
    [M, P] = osid1(test4.A, 25, 4, 1, 1, seed);
    run_osid_test(test4.A, M, P, 28, 1, 0.1);
    [M, P] = osid1(test4.A, 5, 5, 1, 0, seed);
    run_osid_test(test4.A, M, P, 3, 0, 0.1);

    % Wide matrices.
    [M, P] = osid1(test5.A, 27, 3, 1, 1, seed);
    run_osid_test(test5.A, M, P, 29, 1, 0.1);
    [M, P] = osid1(test5.A, 25, 4, 1, 1, seed);
    run_osid_test(test5.A, M, P, 28, 1, 0.1);
    [M, P] = osid1(test5.A, 5, 5, 1, 0, seed);
    run_osid_test(test5.A, M, P, 3, 0, 0.1);

    % osid1
    % Tall matrices.
    [M, P] = osid2(test4.A, 27, 3, 1, 1, seed);
    run_osid_test(test4.A, M, P, 29, 1, 0.1);
    [M, P] = osid2(test4.A, 25, 4, 1, 1, seed);
    run_osid_test(test4.A, M, P, 28, 1, 0.1);
    [M, P] = osid2(test4.A, 5, 5, 1, 0, seed);
    run_osid_test(test4.A, M, P, 3, 0, 0.1);

    % Wide matrices.
    [M, P] = osid2(test5.A, 27, 3, 1, 1, seed);
    run_osid_test(test5.A, M, P, 29, 1, 0.1);
    [M, P] = osid2(test5.A, 25, 4, 1, 1, seed);
    run_osid_test(test5.A, M, P, 28, 1, 0.1);
    [M, P] = osid2(test5.A, 5, 5, 1, 0, seed);
    run_osid_test(test5.A, M, P, 3, 0, 0.1);
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

function[] =  run_osid_test(A, M, P, k, axis, test_tol)

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
        [A_id_ref, ~, ~] = reference_osid(A, k, axis);
        err_ref = norm(A - A_id_ref, 'fro');
        rel_err = (err_rand - err_ref) / norm(A, 'fro');
    end
    disp(rel_err)
    assert(rel_err < test_tol)
end
