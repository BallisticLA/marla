function [Q] = rf1(A, k, p, s)
    %{
    Standard routine for constructing a matrix Q of size 
    (size(A, 2), k) where range(Q) is "reasonably" well aligned with 
    the span of the top k right singular vectors of A.

    Uses power iteration technique (recommended for cases with slow decay 
    of singular values of A), controlled by parameter p. Each power
    iteratin step requires two passes over A.

    s is an int or RandomStream. It controls all random number generation.

    Uses QR decomposition for insuring orthogonality of the columns of
    sketch Q.
    %}
    s = MarlaRandStream(s);
    class_A = class(A);
    [~, n] = size(A);
    % By default, a Gaussian random sketching matrix is used.
    % Alternative choices are present in '../../utils/sketching_operators'.
    Omega = randn(s, n, k, class_A);
    [Q, ~] = qr(A * Omega, 0);

    for j = 1 : p
        [Q, ~] = qr(A' * Q, 0);
        [Q, ~] = qr(A * Q, 0);
    end
end