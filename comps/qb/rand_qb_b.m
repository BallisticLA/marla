function [Q, B, log] = rand_qb_b(A, block_size, tol, k, p, s, logging)
%{
    Iteratively build an approximate QB factorization of A,
    which terminates once one of the following conditions
    is satisfied
        (1)  || A - Q B ||_Fro / ||A|| <= tol
    or
        (2) Q has k columns.
    or
        (2) || A - Q B ||_Fro / ||A|| has increased in comparison to previous
        iteration's result (unlikely case). 

    Each iteration involves sketching A from the right by a sketching
    matrix with "block_size" columns, and using "p" power iterations. 
    
    Parameters
    ----------
    A : matrix
        Data matrix to approximate.

    block_size : int
        The block size in this blocked QB algorithm. Add this many columns
        to Q at each iteration (except possibly the final iteration).

    tol : float
        Terminate if ||A - Q B||_Fro / ||A|| <= tol. Setting k = min(size(A)) is a
        valid way of ensuring ||A - Q B||_Fro / ||A|| <= tol on exit.

    k : int
        Terminate if size(Q, 2) == k. Assuming k < rank(A), setting tol=0 is a
        valid way of ensuring size(Q, 2) == k on exit.

    p : int
        Number of power iterations used. Using this algorithm version
        implies increase of the number of passes over matrix A by 2 with each
        power iteration. 

    s : int or RandomStream
        Controls all random number generation

    logging : struct array 
        Parameter for logging different levels of detailed information.
        Contains two fields:
        Depth - describes on how many subroutines levels should the info
        be logged (0 for none, 1 for just the main routine, etc).
        Span - describes which types of info should be logged (0 for none,
        1 for timings, 2 for timings + errors estimates, 3 for unusual 
        behaviors).
    
    Returns
    -------
    Q : matrix
        Has the same number of rows of A, and orthonormal columns.

    B : matrix
        Has the same number of columns of A.

    log : structure array
        Holds fields with logged information on routine - fields depend on
        subroutines used.
    
    Notes
    -----
    The number of columns in Q increase by "block_size" at each iteration, unless
    that would bring size(Q, 2) > k. In that case, the final iteration only
    adds enough columns to Q so that size(Q, 2) == k.
    We perform p steps of subspace iteration for each
    block of the QB factorization. We stabilize subspace iteration with
    QR factorization at each step.
    References
    ----------
    Algorithm 2 from YGL:2018.
    ----------
    Important note:
    before calling this routine, use:
    addpath('../utils') - for MatrlaRandStream.m
    addpath('../randgefinders') - for different versions of rangefinders.
%}

    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        %disp('Optional parameter for logging detailed information has not been passed.'); 
    else
        log_present = 1;
        logging.depth = logging.depth - 1;
        log.t_sketch = [];
        log.t_power_iter = [];
        log.t_reorthog = [];
        log.relative_error = [];
        log.Bi_norm = [];
    end

    s = MarlaRandStream(s);
    norm_A = norm(A, 'fro');
    if logging.span >= 2 && log_present
        log.A_norm = norm(A, 'fro');
    end
    % Early termination check on an empty input. 
    if norm_A == 0
        fprintf('The input matrix is empty.');
        return
    end
    % Setting initial error to zero.
    approximation_error = 0;
    class_A = class(A);
    [m, n] = size(A);
    norm_B = 0;
    % Pre-initialization of output matrices. 
    Q = zeros(m, 0, class_A);
    B = zeros(0, n, class_A);
    % Iterative stage.
    for i = 1 : ceil((k / block_size))
        if size(B, 1) + block_size > k
            block_size = k - size(B, 1);
        end
        % Constructing a sketch for current iteration. 
        if log_present, tic, end
        Omega_i = randn(s, n, block_size, class_A);
        [Q_i, ~] = qr((A * Omega_i) - (Q * (B * Omega_i)), 0);
        if log_present, log.t_sketch = [log.t_sketch; toc]; end

        if log_present, tic, end
        % Power iterations. 
        for j = 1 : p
            [Q_i, ~] = qr(A' * Q_i - B' * (Q' * Q_i), 0);
            [Q_i, ~] = qr(A * Q_i - Q * (B * Q_i), 0);
        end
        if log_present, log.t_power_iter = [log.t_power_iter; toc]; end

        % Ensuring orthogonalization of Q. 
        if log_present, tic, end
        [Q_i, ~] = qr(Q_i - (Q * (Q' * Q_i)), 0);
        if log_present, log.t_reorthog = [log.t_reorthog; toc]; end
        % Deterministic computation of B_i.
        B_i = Q_i' * A;
        % Approximation error check.
        norm_B = hypot(norm_B, norm(B_i, 'fro'));
        prev_error = approximation_error;
        approximation_error = sqrt(abs(norm_A - norm_B) * (norm_A + norm_B)) / norm_A; 
        % Handling the round-off error accumulation.
        if (i > 1) && (approximation_error > prev_error)
            break
        end
        % Output update. 
        Q = [Q, Q_i]; %#ok<AGROW>
        B = [B; B_i]; %#ok<AGROW>
        
        % Logging approximation error.
        if logging.span >= 2 && log_present
            log.relative_error = [log.relative_error; approximation_error];
            log.Bi_norm = [log.Bi_norm; norm_B];
        end

        % Condition of reaching tolerance. 
        if approximation_error < tol       
            break;
        end
    end
end
