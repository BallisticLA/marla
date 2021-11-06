function [A, S] = gen_test_mat(m, n, k, spectrum)
    Buf = randn(m, k);
    [U, ~] = qr(Buf, 0);
    Buf = randn(n, k);
    [V, ~] = qr(Buf, 0);
    if isscalar(spectrum)
        spectrum = abs(randn(1, k));
        spectrum = sort(spectrum,'descend');
    end
    S = spdiags(spectrum', 0, k, k);
    A = U * S * V';
end