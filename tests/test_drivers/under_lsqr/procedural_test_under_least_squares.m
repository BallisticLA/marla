% Throws an asserion error if something went wrong.
function [] = procedural_test_under_least_squares()
    seed0 = 3893874;
    addpath('../matrix_generators/');

    % linspace_spec
    m = 1000;
    n = 100;
    cond_num = 1e5;
    spectrum = linspace(cond_num^0.5, cond_num^(-0.5), n);
    seed = MarlaRandStream(seed0);
    A = gen_test_mat(m, n, n, spectrum, seed);
    test1.A = A;
    test1.s = spectrum;
    y0 = randn(seed, m, 1);
    test1.c = A' * y0;
    test1.y_opt = A' \ test1.c;
    
    
    % logspace_spec
    m = 1000;
    n = 100;
    cond_num = 1e5;
    spectrum = logspace(log10(cond_num)/2, -log10(cond_num)/2, n);
    seed = MarlaRandStream(seed0);
    A = gen_test_mat(m, n, n, spectrum, seed);
    test2.A = A;
    test2.s = spectrum;
    y0 = randn(seed, m, 1);
    test2.c = A' * y0;
    test2.y_opt = A' \ test2.c;
   
    % lowrank linspace
    m = 1000;
    n = 100;
    cond_num = 1e5;
    rk = 80;
    spectrum = linspace(cond_num^0.5, cond_num^(-0.5), rk);
    seed = MarlaRandStream(seed0);
    A = gen_test_mat(m, n, rk, spectrum, seed);
    test3.A = A;
    test3.s = spectrum;
    y0 = randn(seed, m,1);
    test3.c = A' * y0;
    test3.y_opt = A' \ test3.c;
    
    
    test_spu1(test1, test2, test3);
end


function[] = test_spu1(test1, test2, test3)
    addpath('../../drivers/least_squares/');
    seed = 998765;
        
    % linspace spec
    test1.y_approx = spu1(test1.A, test1.c, 3, 1e-12, 50, 0, seed);
    run_test(test1, 1e-6)

    % logspace spec
    test2.y_approx = spu1(test2.A, test2.c, 3, 1e-12, 50, 0, seed);
    run_test(test2, 1e-6);

    % lowrank linspace spec
    test3.y_approx = spu1(test3.A, test3.c, 3, 1e-12, 50, 0, seed);
    run_test(test3, 1e-6);
end

function[] = run_test(self, test_tol)
    test_objective(self, test_tol);
    test_residual(self, test_tol);
    % test_delta_y(self, test_tol);  % skip this (here, and in Python)
end


function[] = test_delta_y(self, tol)
    % ||y - y_opt|| <= tol
    delta_y = self.y_opt - self.y_approx;
    nrm = norm(delta_y, 2) / (1 + min(norm(self.y_opt, 2), norm(self.y_approx, 2)));
    disp(nrm); 
    assert(nrm <= tol);
end

function[] = test_objective(self, tol)
    % ||y|| <= ||y_opt|| + tol
    nrm_opt = norm(self.y_opt, 2);
    nrm_approx = norm(self.y_approx, 2);
    assert(nrm_approx <= (1 + tol)*nrm_opt);
end

function[] = test_residual(self, tol)
    % ||A' y - c|| <= tol
    res = self.A' * self.y_approx - self.c;
    nrm = norm(res, 2);
    assert(nrm <= tol);
end
