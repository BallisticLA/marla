function [Q] = rangefinder_alter(A, k, p)
    %{
    Routine for constructing a matrix Q of size 
    (size(A, 2), k) where range(Q) is "reasonably" well aligned with 
    the span of the top k right singular vectors of A.

    Uses power iteration technique (recommended for cases with slow decay 
    of singular values of A), controlled by parameter p. Each power
    iteratin step requires two passes over A.

    Uses QR decomposition for insuring orthogonality of the columns of
    sketch Q.

    Differs from standard scheme in using LU factorization for
    stabilization of Q during power iteration steps. 
    %}
    class_A = class(A);
    [~, n] = size(A);
    % By default, a Gaussian random sketching matrix is used.
    % Alternative choices are present in '../../utils/sketching_operators'.
    Omega = randn(n, k, class_A);
    [Q, ~] = qr(A * Omega, 0);

    for j = 1 : p
    
        [Q, ~] = lu(A' * Q);
        
        if j < power 
            [Q, ~] = lu(A * Q);
        else
            [Q, ~] = qr(A * Q, 0);
        end
    end
end