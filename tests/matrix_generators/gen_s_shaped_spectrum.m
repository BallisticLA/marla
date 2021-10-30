function [A, s] = gen_s_shaped_spectrum(m, n, k)
    spectrum = 0.0001+1./(1 + exp((1:k)-30));
    [A, s] = gen_test_mat(m, n, k, spectrum);
end