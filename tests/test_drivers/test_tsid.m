% Riley note: trhese
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

    test_simple_exact(test1, test3);
    test_simple_approx(test4, test5);
end

function[] = test_simple_exact(test1, test3)
    addpath('../../drivers/interpolative/')
    [Z, Is, X, Js] = tsid1(test1.A, test1.k, 0, 1);
    run_tsid_test(test1.A, Z, Is, X, Js, 1e-12);
    
    [Z, Is, X, Js] = tsid1(test1.A, test1.k, 5, 1);
    run_tsid_test(test1.A, Z, Is, X, Js, 1e-12);
    % ^ Once we add more tests (e.g., that A[Is,Js]
    %   is nonsingular), this case might fail.
    
    [Z, Is, X, Js] = tsid1(test3.A, test3.k, 0, 1);
    run_tsid_test(test3.A, Z, Is, X, Js, 1e-12);
end

function[] = test_simple_approx(test4, test5)
    addpath('../../drivers/interpolative/')

    % Tall matrices.
    [Z, Is, X, Js] = tsid1(test4.A, 27, 3, 2);
    run_tsid_test(test4.A, Z, Is, X, Js, 0.05);
    
    [Z, Is, X, Js] = tsid1(test4.A, 25, 4, 2);
    run_tsid_test(test4.A, Z, Is, X, Js, 0.05);

    % Wide matrices.
    [Z, Is, X, Js] = tsid1(test5.A, 27, 3, 2);
    run_tsid_test(test5.A, Z, Is, X, Js, 0.05);
    
    [Z, Is, X, Js] = tsid1(test5.A, 25, 4, 2);
    run_tsid_test(test5.A, Z, Is, X, Js, 0.05);
end

function[] =  run_tsid_test(A, Z, I, X, J, test_tol)
    A_id = Z * (A(I, J) * X);
    err_rand = norm(A - A_id, 'fro');
    if test_tol < 1e-8
        rel_err = err_rand / norm(A, 'fro');
    else
        k = max(size(I));  % should be size(I, 1), but let's be safe.
        [A_id_ref, ~, ~] = reference_osid(A, k, 1);
        err_ref = norm(A - A_id_ref, 'fro');
        rel_err = (err_rand - err_ref) / norm(A, 'fro');
    end
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