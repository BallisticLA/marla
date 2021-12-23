% Throws an asserion error if something went wrong.
function [] = procedural_test_qb()
    seed = 42;
    
    addpath('../matrix_generators/');
    % tall_low_exact_rank
    test1.k = 15;
    test1.A = gen_test_mat(200, 50, test1.k, 0, seed); 

    % wide_low_exact_rank
    test2.k = 15;
    test2.A = gen_test_mat(50, 200, test2.k, 0, seed); 
    
    % tall_full_exact_rank
    test3.k = 50;
    test3.A = gen_test_mat(200, 50, test3.k, 0, seed); 
    
    % wide_full_exact_rank
    test4.k = 50;
    test4.A = gen_test_mat(50, 200, test4.k, 0, seed); 
    
    % fast_decay_low_exact_rank
    test5.k = 15;
    test5.A = gen_exp_spectrum(200, 50, test5.k, 2, seed); 
    
    % slow_decay_low_exact_rank
    test6.k = 15;
    test6.A = gen_exp_spectrum(200, 50, test6.k, 150, seed); 
    
    % s_shaped_decay_low_exact_rank
    test7.k = 15;
    test7.A = gen_s_shaped_spectrum(200, 50, test7.k, seed); 

    % Test rand_QB.
    test_rand_qb(test1, test2, test3, test4, test5, test6, test7);
    % Test rand_QB_B.
    test_rand_qb_b(test1, test2, test3, test4, test5, test6, test7);
    % Test rand_QB_B_PE.
    test_rand_qb_b_pe(test1, test2, test3, test4, test5, test6, test7);
    % Test rand_QB_SP.
    test_rand_qb_sp(test1, test2, test3, test4, test5, test6, test7);
end

function[] = test_rand_qb(test1, test2, test3, test4, test5, test6, test7)
    seed = 24;
    addpath('../../comps/qb/');
    alg_tol = 1e-8;
    p = 1;
    % Run the above algorithm on tall matrices, wide matrices,
    % and matrices with varying types of spectral decay.
    %
    % In all cases we set the target rank of the approximation matrix
    % to the actual rank of the data matrix.
    [Q, B] = rand_qb(test1.A, test1.k, p, seed);
    run_exact(test1.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb(test2.A, test2.k, p, seed);
    run_exact(test2.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb(test5.A, test5.k, p, seed);
    run_exact(test5.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb(test6.A, test6.k, p, seed);
    run_exact(test6.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb(test7.A, test7.k, p, seed);
    run_exact(test7.A, Q, B, alg_tol);

    % Cases with full-rank factorizations.
    [Q, B] = rand_qb(test3.A, test3.k, p, seed);
    run_exact(test3.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb(test4.A, test3.k, p, seed);
    run_exact(test4.A, Q, B, alg_tol);

    % Tests with overestimated rank.
    [Q, B] = rand_qb(test1.A, floor(test1.k * 1.2), p, seed);
    run_estimated(test1.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb(test2.A, floor(test2.k * 1.2), p, seed);
    run_estimated(test2.A, Q, B, alg_tol);
end

function[] = test_rand_qb_b(test1, test2, test3, test4, test5, test6, test7)
    seed = 24;
    addpath('../../comps/qb/');
    alg_tol = 1e-8;
    test_tol = 1e-8;
    block_size = 5;
    p = 1;
    % Run the above algorithm on tall matrices, wide matrices,
    % and matrices with varying types of spectral decay.
    %
    % In all cases we set the target rank of the approximation matrix
    % to the actual rank of the data matrix.
    [Q, B] = rand_qb_b(test1.A, block_size, test_tol, test1.k, p, seed);
    run_exact(test1.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b(test2.A, block_size, test_tol, test2.k, p, seed);
    run_exact(test2.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b(test5.A, block_size, test_tol, test5.k, p, seed);
    run_exact(test5.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b(test6.A, block_size, test_tol, test6.k, p, seed);
    run_exact(test6.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b(test7.A, block_size, test_tol, test7.k, p, seed);
    run_exact(test7.A, Q, B, alg_tol);

    % Cases with full-rank factorizations.
    [Q, B] = rand_qb_b(test3.A, block_size, test_tol, test3.k, p, seed);
    run_exact(test3.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b(test4.A, block_size, test_tol, test4.k, p, seed);
    run_exact(test4.A, Q, B, alg_tol);

    % Tests with overestimated rank.
    [Q, B] = rand_qb_b(test1.A, block_size, test_tol, floor(test1.k * 1.2), p, seed);
    run_estimated(test1.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b(test2.A, block_size, test_tol, floor(test2.k * 1.2), p, seed);
    run_estimated(test2.A, Q, B, alg_tol);    
end

function[] = test_rand_qb_b_pe(test1, test2, test3, test4, test5, test6, test7)
    seed = 24;
    addpath('../../comps/qb/');
    alg_tol = 1e-8;
    test_tol = 1e-8;
    block_size = 5;
    p = 1;
    % Run the above algorithm on tall matrices, wide matrices,
    % and matrices with varying types of spectral decay.
    %
    % In all cases we set the target rank of the approximation matrix
    % to the actual rank of the data matrix.
    [Q, B] = rand_qb_b_pe(test1.A, block_size, test_tol, test1.k, p, seed);
    run_exact(test1.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b_pe(test2.A, block_size, test_tol, test2.k, p, seed);
    run_exact(test2.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b_pe(test5.A, block_size, test_tol, test5.k, p, seed);
    run_exact(test5.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b_pe(test6.A, block_size, test_tol, test6.k, p, seed);
    run_exact(test6.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b_pe(test7.A, block_size, test_tol, test7.k, p, seed);
    run_exact(test7.A, Q, B, alg_tol);

    % Cases with full-rank factorizations.
    [Q, B] = rand_qb_b_pe(test3.A, block_size, test_tol, test3.k, p, seed);
    run_exact(test3.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b_pe(test4.A, block_size, test_tol, test4.k, p, seed);
    run_exact(test4.A, Q, B, alg_tol);

    % Tests with overestimated rank.
    [Q, B] = rand_qb_b_pe(test1.A, block_size, test_tol, floor(test1.k * 1.2), p, seed);
    run_estimated(test1.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_b_pe(test2.A, block_size, test_tol, floor(test2.k * 1.2), p, seed);
    run_estimated(test2.A, Q, B, alg_tol);  
end

function[] = test_rand_qb_sp(test1, test2, test3, test4, test5, test6, test7)
    seed = 24;
    addpath('../../comps/qb/');
    alg_tol = 1e-8;
    % Run the above algorithm on tall matrices, wide matrices,
    % and matrices with varying types of spectral decay.
    %
    % In all cases we set the target rank of the approximation matrix
    % to the actual rank of the data matrix.
    [Q, B] = rand_qb_sp(test1.A, test1.k, seed);
    run_exact(test1.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_sp(test2.A, test2.k, seed);
    run_exact(test2.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_sp(test5.A, test5.k, seed);
    run_exact(test5.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_sp(test6.A, test6.k, seed);
    run_exact(test6.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_sp(test7.A, test7.k, seed);
    run_exact(test7.A, Q, B, alg_tol);

    % Cases with full-rank factorizations.
    [Q, B] = rand_qb_sp(test3.A, test3.k, seed);
    run_exact(test3.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_sp(test4.A, test4.k, seed);
    run_exact(test4.A, Q, B, alg_tol);

    % Tests with overestimated rank.
    [Q, B] = rand_qb_sp(test1.A, floor(test1.k * 1.2), seed);
    run_estimated(test1.A, Q, B, alg_tol);
    
    [Q, B] = rand_qb_sp(test2.A, floor(test2.k * 1.2), seed);
    run_estimated(test2.A, Q, B, alg_tol);
end
% 

function[] = run_exact(A, Q, B, tol)
    test_exact(A, Q, B, tol);
    test_valid_onb(Q, tol);
    test_exact_B(A, Q, B, tol);
    test_exact_B(A, Q, B, tol);
    test_exact_rank_B(B);
end

% Same as above - rank(B) test.
function[] = run_estimated(A, Q, B, tol)
    test_exact(A, Q, B, tol);
    test_valid_onb(Q, tol);
    test_exact_B(A, Q, B, tol);
    test_exact_B(A, Q, B, tol);
end

function[] = test_exact(A, Q, B, tol)
    delta = A - Q * B;
    nrm = norm(delta, 'fro');
    assert(nrm <= tol);
end

function[] = test_valid_onb(Q, tol)
    gram = Q' * Q;
    delta = gram - eye(size(gram));
    nrm = norm(delta, 'fro');
    assert(nrm <= tol);
end

function[] = test_exact_B(A, Q, B, tol)
    delta = B - Q' * A;
    nrm = norm(delta, 'fro');
    assert(nrm <= tol);
end

function[] = test_exact_rank_B(B)
   assert(rank(B) == min(size(B))); 
end