function [Q, log] = rf1(A, k, p, seed, logging)
    %{
    Standard routine for constructing an orthonormal sketch Q of size 
    (size(A, 2), k) where range(Q) is "reasonably" well aligned with 
    the span of the top k right singular vectors of A.

    s is an int or RandomStream. It controls all random number generation.

    Uses QR decomposition for insuring orthogonality of the columns of
    sketch Q.
    %}
    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        log.status = 'Optional parameter for logging detailed information has not been passed.'; 
    else
        log_present = 1;
    end

    % Set to a default value. 
    passes_per_stab = 1;

    seed = MarlaRandStream(seed);
    [Omega, log] = rs1(A, k, p, passes_per_stab, seed, logging);
    [Q, ~] = qr(A * Omega, 0);
end