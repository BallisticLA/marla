function [res, log] = spo3(A, b, sampling_factor, tol, iter_lim, use_chol, logging, s)
    %{
    A sketch-and-precondition approach to overdetermined ordinary least
    squares. This implementation uses QR (when use_chol=false) or Cholesky
    (when use_chol=true) to obtain the preconditioner. This implementation
    uses LSQR for the iterative method.

    Before starting LSQR, we run a basic sketch-and-solve (essentially for
    free) to obtain a solution x_ske. If ||A * x_ske - b||_2 < ||b||_2,
    then we initialize LSQR at x_ske.
 
    This implementation assumes A is full rank.

    References
    ----------
    This implementation was inspired by Blendenpik (AMT:2010). The differences
    relative to the official Blendenpik algorithm [AMT:2010, Algorithm 1] are

        (1) We make no assumption on the distribution of the sketching matrix
            which may be to form the preconditioner. Blendenpik only used
            SRTTs (Walsh-Hadamard, discrete cosine, discrete Hartley).
            However, this specific implementation uses random Gaussian
            sketching matrix. 

        (2) We let the user specify the exact embedding dimension, as
            floor(self.oversampling_factor * A.shape[1]).

        (3) We do not zero-pad A with additional rows. Such zero padding
            might be desirable to facilitate rapid application of SRTT
            sketching operators. It is possible to implement an SRTT operator
            so that it performs zero-padding internally.

        (4) We do not perform any checks on the quality of the preconditioner.

        (5) We initialize the iterative solver (LSQR) at the better of the two
            solutions given by either the zero vector or the output of
            sketch-and-solve.

    Note:
    This implementation has the option of logging detailed information
    on runtime and the rate at which (preconditioned) normal equation
    error decays while LSQR runs. Controlled by passing a boolean parameter
    "logging".

    Important note:
    Before running, use: 
    addpath('../../Utils/Sketching_Operators')
    %}
    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        log.status = 'Optional parameter for logging detailed information has not been passed.'; 
    else
        log_present = 1;
        logging.depth = logging.depth - 1;
    end

    s = MarlaRandStream(s);
    
    [num_rows, num_cols] = size(A);
    d = lstsq_dim_checks(sampling_factor, num_rows, num_cols);
    
    % Sketch the data matrix.
    if log_present, tic, end
    % By default, a SJLT sketching matrix is used.
    % Alternative choices are present in '../../Utils/Sketching_Operators'.
    Omega = sjlt(double(d), double(num_rows), double(8), s);
    A_ske = Omega * A;
    if log_present, log.t_sketch = toc; end
    
    % Factor the sketch.
    if log_present, tic, end 
    if ~use_chol
        [Q, R] = qr(A_ske, 0);
    else
        % It would be better to call "syrk" from BLAS when forming the
        % Gram matrix
        R = chol(A_ske' * A_ske, 'upper');
    end
    if log_present, log.t_factor = toc; end
    
    % Sketch-and-solve type preprocessing.
    if log_present, tic, end
    b_ske = Omega * b;
    if ~use_chol
        z_ske = Q' * b_ske;
    else
        z_ske = R' \ (A_ske' * b_ske);
    end
    x_ske = R \ z_ske;
    if norm(A * x_ske - b, 'fro') >= norm(b)
        z_ske = zeros(num_cols, 1);
    end
    if log_present, log.t_presolve = toc; end
    
    % Iterative phase.
    if log_present, tic, end
    [res, ~, ~, iters, resvec, ~] = lsqr(A, b, tol, iter_lim, R, eye(size(R, 1)), x_ske);
    if log_present, log.t_iterate = toc; end    
    
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
        ar0 = A' * b;
        ar0 = R \ ar0;
        ar0norm = norm(ar0, 'fro');
        log.x = res;
        log.errors = [ar0norm; resvec];
    end
end
