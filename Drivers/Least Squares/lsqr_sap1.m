function [res, log] = lsqr_sap1(A, b, sampling_factor, tol, iter_lim, logging)
    %{
    A sketch-and-precondition approach to overdetermined ordinary least
    squares. This implementation uses QR to obtain the preconditioner and
    it uses LSQR for the iterative method.

    Before starting LSQR, we run a basic sketch-and-solve (for free, given
    our QR decomposition of the sketched data matrix) to obtain a solution
    x_ske. If ||A * x_ske - b||_2 < ||b||_2, then we initialize LSQR at x_ske.
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
    %}
    if ~logging
        disp('Optional parameter for logging detailed information has not been passed.'); 
    end
    
    [num_rows, num_cols] = size(A);
    d = dim_checks(sampling_factor, num_rows, num_cols);
    
    % Sketch the data matrix.
    if logging, tic, end
    % By default, a SJLT sketching matrix is used.
    % Alternative choices are present in '../../Utils/Sketching_Operators'.
    addpath('../../Utils/Sketching_Operators')
    Omega = sjlt(d, num_rows, 8);
    A_ske = Omega * A;
    if logging, log.t_sketch = toc; end
    
    % Factor the sketch.
    if logging, tic, end 
    [Q, R] = qr(A_ske, 0);
    if logging, log.t_factor = toc; end
    
    % Sketch-and-solve type preprocessing.
    if logging, tic, end
    b_ske = Omega * b;
    x_ske = R \ (Q' * b_ske);
    disp(x_ske);
    x0 = int16.empty;
    
    if norm(A * x_ske - b, 'fro') < norm(b)
        x0 = x_ske;
    end
    if logging, log.t_presolve = toc; end
    
    % Iterative phase.
    if logging, tic, end
    [res, ~, ~, iters, resvec, ~] = lsqr(A, b, tol, iter_lim, R, eye(size(R, 1)), x0);
    if logging, log.t_iterate = toc; end    
    
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
        ar0 = A' * b;
        ar0 = R \ ar0;
        ar0norm = norm(ar0, 'fro');
        log.x = res;
        log.arnorms = [ar0norm; resvec];
    end
end

% Helper routine. 
function [d] = dim_checks(sampling_factor, num_rows, num_cols)
    assert(num_rows >= num_cols);
    d = cast((sampling_factor * num_cols), 'uint8');
    if d > num_rows
        fprintf(['The embedding dimension "d" should not be larger than the', ... 
        'number of rows of the data matrix. Here, an embedding dimension', ...
        'of d=%d has been requested for a matrix with only %d rows.', ...
        'We will proceed by setting d=%d. This parameter choice will', ...
        'result in a very inefficient algorithm!'], d, num_rows, num_rows);
        d = num_rows; 
    end
    assert(d >= num_cols);
end