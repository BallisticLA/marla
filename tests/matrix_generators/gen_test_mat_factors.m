function [U, spectrum, Vt] = gen_test_mat_factors(m, n, k, spectrum, s)
    addpath('../../../utils/');
    s = MarlaRandStream(s);
    Buf = randn(s, m, k);
    [U, ~] = qr(Buf, 0);
    Buf = randn(s, n, k);
    [V, ~] = qr(Buf, 0);
    Vt = V';
    if isscalar(spectrum)
        spectrum = abs(randn(s, 1, k));
        spectrum = sort(spectrum,'descend');
    end
end