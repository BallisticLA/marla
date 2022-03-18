function [Omega, log] = rs1(A, k, p, passes_per_stab, seed, logging)
    %{
    Routine for constructiong a random sketching operator Omega, using an
    optimized power iteration scheme that allows for any nuber of passes
    over initial data matrix A. 

    Power iteration technique (recommended for cases with slow decay 
    of singular values of A) is controlled by parameter 'p' (>= 0). Additionally,
    stabilization (used to avoid accumulation of numerical error) is
    controlled by parameter 'passes_per_stab' (>= 1). By default, utilizes
    LU for stabilization (alternatively, can use TSQR).

    seed is an int or RandomStream. It controls all random number generation.
    %}
    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        log.status = 'Optional parameter for logging detailed information has not been passed.'; 
    else
        log_present = 1;
    end

    class_A = class(A);
    [m, n] = size(A);
    p_done = 0;

    % By default, a Gaussian random sketching matrix is used.
    % Alternative choices are present in
    % '../../utils/sketching_operators'.
    if log_present, tic, end
    if mod(p, 2) == 0
        Omega = randn(seed, n, k, class_A);     
    else
        Omega = A' * randn(seed, m, k, class_A);
        p_done = p_done + 1;
        if mod(p_done, passes_per_stab) == 0
            [Omega, ~] = lu(Omega);
        end
    end
    
    while (p - p_done) >= 2
        Omega = A * Omega;
        p_done = p_done + 1;
        if mod(p_done, passes_per_stab) == 0
            [Omega, ~] = lu(Omega);
        end
        Omega = A' * Omega;
        p_done = p_done + 1;
        if mod(p_done, passes_per_stab) == 0
            [Omega, ~] = lu(Omega);
        end
    end
    if log_present, log.t_sketch = toc; end
end
