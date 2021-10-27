function [Q] = skecth_construction_PE(A, k, p)
    %{
    Pass-Efficient routine for constructing a matrix Q of size 
    (size(A, 2), k) where range(Q) is "reasonably" well aligned with 
    the span of the top k right singular vectors of A.

    Uses power iteration technique (recommended for cases with slow decay 
    of singular values of A), controlled by parameter p.

    Uses QR decomposition for insuring orthogonality of the columns of
    sketch Q.

    Allows for any number of passes over the initial matrix A.
    Uses LU factorization for stabilization during power iterations step. 
    %}
    class_A = class(A);
    [m, n] = size(A);
    v = 2 * p + 1;

    % Odd number of passes over A.
    if(mod(v, 2) == 0)
        % By default, a Gaussian random sketching matrix is used.
        % Alternative choices are present in '../Sketching_Operators'
        Omega = randn(m, k, class_A);
        
        if (v > 2)
            [Q, ~] = lu(A' * Omega);
        else
            [Q, ~] = qr(A' * Omega, 0);
        end
    % Even number of passes over A.
    else
        % By default, a Gaussian random sketching matrix is used.
        % Alternative choices are present in '../Sketching_Operators'
        Q = randn(n, k, class_A);
        if p == 0
            [Q, ~] = qr(A * Q, 0);
        end
    end

    for i = 1 : p

        [Q, ~] = lu(A * Q);
        if i == p
            [Q, ~] = qr(A' * Q, 0);
        else
            [Q, ~] = lu(A' * Q);
        end
    end
end