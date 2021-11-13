%{
    NOTE:
    Set relative error instead of absolute error due to the way rand_QB
    operates. 

    Occasionally, test_valid_Onv fails for rand_qb_b 
    with true error ~2.1237e-08 vs error tol 1.0e-8. 
    -on test_fr

    Occasionally, test_abs_fro_error fails for rand_QB_B 
    with true error ~3.2536e-08 vs error tol 1.0000e-08.
    -On test_fr
    
    To investigate.
%}
function [] = test_svd()
    addpath('../matrix_generators/');
    seed = 59;

    % tall_low_exact_rank
    test1.k = 15;
    test1.A = gen_test_mat(200, 50, test1.k, 0, seed); 
    test1.S_exact = svd(test1.A); 

    % wide_low_exact_rank
    test2.k = 15;
    test2.A = gen_test_mat(50, 200, test2.k, 0, seed); 
    test2.S_exact = svd(test2.A); 

    % tall_full_exact_rank
    test3.k = 50;
    test3.A = gen_test_mat(200, 50, test3.k, 0, seed); 
    test3.S_exact = svd(test3.A); 

    % wide_full_exact_rank
    test4.k = 50;
    test4.A = gen_test_mat(50, 200, test4.k, 0, seed); 
    test4.S_exact = svd(test4.A); 

    test_fr(test1, test2);
    test_fp_inexact(test3, test4);
    test_fp_exact(test1, test2, test3);
end

function[] = test_fr(test1, test2)
    addpath('../../drivers/');
    %{     
    For a wide matrix and a tall matrix:
    Fixed rank QB algorithm
    Three cases:
    Target rank < exact rank (no oversampling).
    Target rank + oversampling < exact rank
    Target rank + oversampling > exact rank
    %}   
    alg_tol = 1e-8;
    seed = 99;

    % Tall matrix;
    rank = size(test1.S_exact, 1);
    [test1.U, test1.S, test1.V] = svd1(test1.A, rank - 10, 0, 1, alg_tol, 5, seed);
    run_batch(test1, alg_tol, 1e-8);
    [test1.U, test1.S, test1.V] = svd1(test1.A, rank - 10, 5, 1, alg_tol, 5, seed);
    run_batch(test1, alg_tol, 1e-8);
    [test1.U, test1.S, test1.V] = svd1(test1.A, rank - 3, 5, 1, alg_tol, 5, seed);
    run_batch(test1, alg_tol, 1e-8);

    % Wide matrix
    rank = size(test2.S_exact, 1);
    [test2.U, test2.S, test2.V] = svd1(test2.A, rank - 10, 0, 1, alg_tol, 5, seed);
    run_batch(test2, alg_tol, 1e-8);
    [test2.U, test2.S, test2.V] = svd1(test2.A, rank - 10, 5, 1, alg_tol, 5, seed);
    run_batch(test2, alg_tol, 1e-8);
    [test2.U, test2.S, test2.V] = svd1(test2.A, rank - 3, 5, 1, alg_tol, 5, seed);
    run_batch(test2, alg_tol, 1e-8);

end

function[] = test_fp_inexact(test3, test4)
    addpath('../../drivers/');
    seed = 11111;
    % set the relative error tolerance to 0.001
    rank = min(size(test3.A));
    rel_err = 0.001;
    [test3.U, test3.S, test3.V] = svd1(test3.A, rank, 0, 1, rel_err, 5, seed);

    run_batch(test3, rel_err, 1e-8);
    [test3.U, test3.S, test3.V] = svd1(test3.A, rank, 2, 1, rel_err, 5, seed);
    run_batch(test3, rel_err, 1e-8);

    % Wide matrix
    rank = min(size(test4.A));
    [test4.U, test4.S, test4.V] = svd1(test4.A, rank, 0, 1, rel_err, 5, seed);
    run_batch(test4, rel_err, 1e-8);
    [test4.U, test4.S, test4.V] = svd1(test4.A, rank, 2, 1, rel_err, 5, seed);
    run_batch(test4, rel_err, 1e-8);
end

function[] = test_fp_exact(test1, test2, test3)
    addpath('../../drivers/');
    % Set the (relative) target tolerance to 1e-12.
    % Tall matrix (low exact rank)
    seed = 1512;
    rank = min(size(test1.A));
    rel_err = 1e-12;
    [test1.U, test1.S, test1.V] = svd1(test1.A, rank, 0, 1, rel_err, 5, seed);
    run_batch(test1, rel_err, 1e-8);
    [test1.U, test1.S, test1.V] = svd1(test1.A, rank, 2, 1, rel_err, 5, seed);
    run_batch(test1, rel_err, 1e-8);

    % Tall matrix (full rank)
    rank = min(size(test3.A));
    [test3.U, test3.S, test3.V] = svd1(test3.A, rank, 0, 1, rel_err, 5, seed);
    run_batch(test3, rel_err, 1e-8);
    [test3.U, test3.S, test3.V] = svd1(test3.A, rank, 1, 1, rel_err, 5, seed);
    run_batch(test3, rel_err, 1e-8);
    [test3.U, test3.S, test3.V] = svd1(test3.A, rank, 2, 1, rel_err, 5, seed);
    run_batch(test3, rel_err, 1e-8);

    % Wide matrix
    rank = min(size(test2.A));
    [test2.U, test2.S, test2.V] = svd1(test2.A, rank, 0, 1, rel_err, 5, seed);
    run_batch(test2, rel_err, 1e-8);
    [test2.U, test2.S, test2.V] = svd1(test2.A, rank, 2, 1, rel_err, 5, seed);
    run_batch(test2, rel_err, 1e-8);
end


function[] = run_batch(self, test_tol, target_tol)
    test_conformable(self);
    test_valid_onb(self, target_tol);
    test_valid_singvals(self, test_tol);
    if ~(isnan(target_tol))
        test_abs_fro_error(self, target_tol);
    end
end

function[] = test_conformable(self)
    assert(size(self.S, 1) == size(self.U, 2));
    assert(size(self.S, 2) == size(self.V, 2));
    assert(size(self.U, 1) == size(self.A, 1));
    assert(size(self.V, 1) == size(self.A, 2));
end

function[] = test_valid_onb(self, tol)
    U = self.U;
    Vt = self.V';

    gram = U' * U;
    delta = gram - eye(size(gram));
    nrm = norm(delta, 'fro');
    assert(nrm <= tol);

    gram = Vt * Vt';
    delta = gram - eye(size(gram));
    nrm = norm(delta, 'fro');
    assert(nrm <= tol);
end

function[]  = test_valid_singvals(self, test_tol)
    self.S = self.S;
    self.S_exact = self.S_exact;

    assert(size(self.S, 1) <= size(self.S_exact, 1));
    assert(min(size(self.S)) >= 0);
end

function[] = test_abs_fro_error(self, rel_tol)
    A = self.A;
    U = self.U;
    S = self.S;
    Vt = self.V';
    % TODO: change this to relative tolerance
    delta = A - (U * S) * Vt;
    nrm = norm(delta, 'fro');
    % abs_tol = rel_tol * norm(S, 2)
    % ^ Scale by  Frobenius norm of A.
    % assert(nrm <= rel_tol)
    assert(nrm <= rel_tol);
end

