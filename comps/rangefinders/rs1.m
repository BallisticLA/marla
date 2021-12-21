function [Q, log] = rs1(A, k, p, s, logging)
    %{
    Pass-Efficient routine for constructing a matrix Q of size 
    (size(A, 2), k) where range(Q) is "reasonably" well aligned with 
    the span of the top k right singular vectors of A.

    Uses power iteration technique (recommended for cases with slow decay 
    of singular values of A), controlled by parameter p.

    s is an int or RandomStream. It controls all random number generation.

    Uses QR decomposition for insuring orthogonality of the columns of
    sketch Q.

    Allows for any number of passes over the initial matrix A.
    Uses LU factorization for stabilization during power iterations step. 
    %}
    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        log.status = 'Optional parameter for logging detailed information has not been passed.'; 
    else
        log_present = 1;
    end

    class_A = class(A);
    [m, n] = size(A);
    v = 2 * p + 1;
    
    if log_present, tic, end
    % Odd number of passes over A.
    if(mod(v, 2) == 0)
        % By default, a Gaussian random sketching matrix is used.
        % Alternative choices are present in '../Sketching_Operators'
        Omega = randn(s, m, k, class_A);
        
        if (v > 2)
            [Q, ~] = lu(A' * Omega);
        else
            [Q, ~] = qr(A' * Omega, 0);
        end
    % Even number of passes over A.
    else
        % By default, a Gaussian random sketching matrix is used.
        % Alternative choices are present in
        % '../../utils/sketching_operators'.
        Q = randn(s, n, k, class_A);
        if p == 0
            [Q, ~] = qr(A * Q, 0);
        end
    end
    if log_present, log.t_sketch = toc; end

    if log_present, tic, end
    for i = 1 : p
        [Q, ~] = lu(A * Q);
        if i == p
            [Q, ~] = qr(A' * Q, 0);
        else
            [Q, ~] = lu(A' * Q);
        end
    end
    if log_present, log.t_power_iter = toc; end
end