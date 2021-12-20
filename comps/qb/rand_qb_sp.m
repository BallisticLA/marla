function [Q, B, log] = rand_qb_sp(A, k, s, logging)
%{
    "Single-pass" version of QB algorithm.
    Return matrices (Q, B) from a rank-k QB factorization of A.
    
    Parameters
    ----------
    A : matrix
        Data matrix to approximate.

    k : int
        Target rank for the approximation of A: 0 < k < min(size(A)).
        This parameter includes any oversampling. For example, if you
        want to be near the optimal (Eckhart-Young) error for a rank 20
        approximation of A, then you might want to set k=25.

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
        Has size (size(A, 1), k). Columns are orthonormal.

    B : matrix
        Has shape (k, size(A, 2)).

    log : structure array
        Holds fields with logged information on routine - fields depend on
        subroutines used.
    
    ----------
    Avoids using power iterations by definition. 
    Requires computing two sketches - in the provided implementation,
    this stage is not done in a single pass, although it can be.
    References
    ----------
    Initial idea stated in section 5.5 of  Halko, Nathan, Per-Gunnar
    Martinsson, and Joel A. Tropp's paper. Further developed in YGL:2018.
    ----------
    Important note:
    before calling this routine, use:
    addpath('../utils') - for MatrlaRandStream.m
%}
    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        log.status = 'Optional parameter for logging detailed information has not been passed.'; 
    else
        log_present = 1;
        logging.depth = logging.depth - 1;
    end

    s = MarlaRandStream(s);
    class_A = class(A);
    [m, n] = size(A);
    % Sketch construction stage - two sketches are constructed here. 
    if log_present, tic, end
    Omega = randn(s, n, k, class_A);
    Omega_ = randn(s, m, k, class_A);
    Y = A * Omega;
    Y_ = A' * Omega_;
    [Q, ~] = qr(Y, 0);
    [Q_, ~] = qr(Y_, 0);
    if log_present, log.t_sketch = toc; end
    % Stage of computing B deterministically.
    B = (Omega_' * Q) \ (Y_' * Q_) * Q_';

    % Relative error computation
    if logging.span >= 2 && log_present
        log.A_fro_norm = norm(A, 'fro');
        log.absolute_error = norm(A - Q * B, 'fro');
    end
end
