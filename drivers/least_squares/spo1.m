function [x_star, log] = spo1(A, b, sampling_factor, tol, iter_lim, smart_init, logging) 
    %{
    A sketch-and-precondition approach to overdetermined ordinary least
    squares. This implementation uses the SVD to obtain the preconditioner
    and it uses LSQR for the iterative method.

    Before starting LSQR, we run a basic sketch-and-solve (for free, given
    our SVD of the sketched data matrix) to obtain a solution x_ske.
    If ||A x_ske - b||_2 < ||b||_2, then we initialize LSQR at x_ske.
    This implementation does not require that A is full-rank.
    References
    ----------
    This implementation was inspired by LSRN. The differences relative to the
    official LSRN algorithm [MSM:2014, Algorithm 1] are
        (1) We make no assumption on the potential distribution of the sketching operator.
            However, this specific implementation uses random Gaussian
            sketching matrix. 

        (2) We provide the option of intelligently initializing the iterative
            solver (LSQR) with the better of the two solutions given by the
            zero vector and the result of sketch-and-solve.
    
    Note:
    This implementation has the option of logging detailed information
    on runtime and the rate at which (preconditioned) normal equation
    error decays while LSQR runs. Controlled by passing a boolean parameter
    "logging".
    %}
    if logging == 0
        %disp('Optional parameter for logging detailed information has not been passed.'); 
    end

    [num_rows, num_cols] = size(A);
    d = lstsq_dim_checks(sampling_factor, num_rows, num_cols);
    
    % Sketch the data matrix.
    if logging, tic, end
    % By default, a SJLT sketching matrix is used.
    % Alternative choices are present in '../../Utils/Sketching_Operators'.
    addpath('../../Utils/Sketching_Operators')
    Omega = sjlt(d, num_rows, 8);
    A_ske = Omega * A;
    if logging, log.t_sketch = toc; end
    
    % Factor the sketch. 
    %  We also measure the time to scale the right singular vectors
    %   as needed for the preconditioner.
    if logging, tic, end
    [U, S, V] = svd(A_ske, 'econ');
    S = diag(S);
    rank = nnz(S > S(1) * num_cols * eps('double'));
    N = V(:, 1:rank) ./ (S(1:rank)');
    if logging, log.t_factor = toc; end
    
    if smart_init
        % Sketch-and-solve type preprocessing.
        %   This isn't necessarily preferable, because it changes the
        %   norm of b, which affects termination criteria.  
        if logging, tic, end
        b_ske = Omega * b;
        z_ske = (U(:, 1:rank)' * b_ske);
        x_ske = N * z_ske;
        b_remainder = b - A * x_ske;
        if norm(b_remainder, 2) >= norm(b, 2)
            z_ske = zeros(size(A, 2), 1);
        end
        if logging, log.t_presolve = toc; end
        
        % Iterative phase.
        if logging, tic, end
        [x_star, iters, resvec] = cgls(A, b, tol, iter_lim, N, x_ske);
        if logging, log.t_iterate = toc; end 
    else
        % No presolve
       if logging, tic, end
       [x_star, iters, resvec] = cgls(A, b, tol, iter_lim, N, zeros(size(A, 2), 1));
       if logging, log.t_iterate = toc; end 
    end
    
    if logging
        % Record a vector of cumulative times to (1) sketch and factor, and
        % (2) take an individual step in LSQR (amortized!).
        %
        % Amortizing the time taken by a single step of LSQR is reasonable,
        % because convergence behavior can be seen by how the normal
        % equation error decays from one iteration to the next.
        t_setup = log.t_sketch + log.t_factor;
        amortized = linspace(0, log.t_iterate, iters);
        cumulative = t_setup + log.t_presolve + amortized;
        log.times = [t_setup cumulative];
        
        % Record a vector of (preconditioned) normal equation errors. Treat
        % the zero vector as a theoretically valid initialization point which
        % we would use before the "solve" phase of "sketch-and-solve".
        ar0 = N' * (A' * b);
        ar0norm = norm(ar0, 'fro');
        log.x = x_star;
        log.errors = [ar0norm; resvec];
    end
end
