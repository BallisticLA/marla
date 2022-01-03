function [x_star, log] = spo1(A, b, sampling_factor, tol, iter_lim, smart_init, logging, seed) 
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

    Input
    -----
    A : matrix
        Data matrix

    b : vector
        Data array

    sampling_factor : float
        Defines the leading dimension of a sketching operator.

    over : int
        Oversampling parameter

    tol : float
        The target error used by the randomized part of the algorithm.

    iter_lim : int
        Maximum number of iterations of lsqr subroutine. 

    smart_init : bool
        Defines initialization of sketch and solve-type preprocessing.

    logging : struct array 
        Parameter for logging different levels of detailed information.
        Contains two fields:
        Depth - describes on how many subroutines levels should the info
        be logged (0 for none, 1 for just the main routine, etc).
        Span - describes which types of info should be logged (0 for none,
        1 for timings, 2 for timings + errors estimates, 3 for unusual 
        behaviors).

    seed: int or RandStream
         Seed for rand stream.

     Output
     ------
     x_star: vector
        Approximation to the solution of the given system.

     log : structure array
         Holds fields with logged information on routine - fields depend on
         subroutines used.

    Important note: 
    ---------------
    Before calling this routine, use:
    addpath('../../utils') - for MatrlaRandStream.m
    addpath('../../Utils/Sketching_Operators')
    %}
    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        log.status = 'Optional parameter for logging detailed information has not been passed.'; 
    else
        log_present = 1;
        logging.depth = logging.depth - 1;
    end

    seed = MarlaRandStream(seed);

    [num_rows, num_cols] = size(A);
    d = lstsq_dim_checks(sampling_factor, num_rows, num_cols);
    
    % Sketch the data matrix.
    if log_present, tic, end
    % By default, a SJLT sketching matrix is used.
    % Alternative choices are present in '../../Utils/Sketching_Operators'.
    Omega = sjlt(d, num_rows, 8, seed);
    A_ske = Omega * A;
    if log_present, log.t_sketch = toc; end
    
    % Factor the sketch. 
    %  We also measure the time to scale the right singular vectors
    %   as needed for the preconditioner.
    if log_present, tic, end
    [U, S, V] = svd(A_ske, 'econ');
    S = diag(S);
    rank = nnz(S > S(1) * num_cols * eps('double'));
    N = V(:, 1:rank) ./ (S(1:rank)');
    if log_present, log.t_factor = toc; end
    
    if smart_init
        % Sketch-and-solve type preprocessing.
        %   This isn't necessarily preferable, because it changes the
        %   norm of b, which affects termination criteria.  
        if log_present, tic, end
        b_ske = Omega * b;
        z_ske = (U(:, 1:rank)' * b_ske);
        x_ske = N * z_ske;
        b_remainder = b - A * x_ske;
        if norm(b_remainder, 2) >= norm(b, 2)
            z_ske = zeros(size(A, 2), 1);
        end
        if log_present, log.t_presolve = toc; end
        
        % Iterative phase.
        if log_present, tic, end
        [x_star, iters, resvec] = cgls(A, b, tol, iter_lim, N, x_ske);
        if log_present, log.t_iterate = toc; end 
    else
        % No presolve
       if log_present, tic, end
       [x_star, iters, resvec] = cgls(A, b, tol, iter_lim, N, zeros(size(A, 2), 1));
       if log_present, log.t_iterate = toc; end 
    end
    
    if log_present
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
