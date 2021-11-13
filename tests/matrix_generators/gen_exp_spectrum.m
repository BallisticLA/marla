function [A, s] = gen_exp_spectrum(m, n, k, t, s)
    s = MarlaRandStream(s);
    spectrum = exp((1 : k) / -t);
    [A, s] = gen_test_mat(m, n, k, spectrum, s);
end