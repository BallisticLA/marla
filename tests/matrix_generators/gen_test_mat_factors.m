function [U, spectrum, V] = gen_test_mat_factors(m, n, k, spectrum, s)
    s = MarlaRandStream(s);
    Buf = randn(s, m, k);
    [U, ~] = qr(Buf, 0);
    Buf = randn(s, n, k);
    [V, ~] = qr(Buf, 0);
    if isscalar(spectrum)
        spectrum = abs(randn(s, 1, k));
        spectrum = sort(spectrum,'descend');
    end
end