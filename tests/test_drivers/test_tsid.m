function[] = test_tsid()
    addpath('../matrix_generators/');

    test1.k = 290;
    test1.A = gen_test_mat(1000, 300, test1.k, 0);  

    test3.k = 30;
    test3.A = gen_test_mat(1000, 300, test3.k, 0); 

    test4.k = 30;
    test4.A = gen_test_mat(100, 30, test4.k, 0); 

    test5.k = 30;
    test5.A = gen_test_mat(30, 100, test5.k, 0); 

    test_simple_exact(test1, test2, test3);
    test_simple_approx(test4, test5);
end

function[] = test_simple_exact(test1, test2, test3)
    addpath('../../drivers/interpolative/')
    [Z, Is, X, Js] = tsid1(test1.A, 290, 0, 1);
    run_tsid_test(test1.A, Z, Is, X, Js, 29, 1e-12);
    [Z, Is, X, Js] = tsid1(test1.A, 290, 5, 1);
    run_tsid_test(test2.A, Z, Is, X, Js, 28, 1e-12);
    [Z, Is, X, Js] = osid1(test3.A, 30, 0, 1);
    run_tsid_test(test3.A, Z, Is, X, Js, 30, 1e-12);
end

function[] = test_simple_approx(test4, test5)
    addpath('../../drivers/interpolative/')

    % Tall matrices.
    [Z, Is, X, Js] = osid1(test4.A, 27, 3, 1, 1);
    run_tsid_test(test4.A, Z, Is, X, Js, 29, 0.05);
    [Z, Is, X, Js] = osid1(test4.A, 25, 4, 1, 1);
    run_tsid_test(test4.A, Z, Is, X, Js, 28, 0.05);

    % Wide matrices.
    [Z, Is, X, Js] = osid1(test5.A, 27, 3, 1, 1);
    run_tsid_test(test5.A, Z, Is, X, Js, 29, 0.05);
    [Z, Is, X, Js] = osid1(test5.A, 25, 4, 1, 1);
    run_tsid_test(test5.A, Z, Is, X, Js, 28, 0.05);
end

function[] =  run_tsid_test(A, Z, I, X, J, k, test_tol)
    A_id = Z * (A(I, J) * X);
    err = norm(A - A_id, 'fro') / norm(A, 'fro');
    assert(err < test_tol);
end
