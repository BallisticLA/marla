function [U, spectrum, V] = gen_test_mat_factors(m, n, k, spectrum)
    Buf = randn(m, k);
    [U, ~] = qr(Buf, 0);
    Buf = randn(n, k);
    [V, ~] = qr(Buf, 0);
    if isscalar(spectrum)
        spectrum = abs(randn(1, k));
        spectrum = sort(spectrum,'descend');
    end
end