function [Q, B, log] = rand_qb(A, k, p, s, logging)
%{
    Return matrices (Q, B) from a rank-k QB factorization of A.
    ----------
    A : matrix
        Data matrix to approximate.
    
    k : int
        Target rank for the approximation of A: 0 < k < min(size(A)).
        This parameter includes any oversampling. For example, if you
        want to be near the optimal (Eckhart-Young) error for a rank 20
        approximation of A, then you might want to set k=25.
    
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
        Has size (size(A, 1), k). Columns are orthonormal.
    
    B : matrix
        Has shape (k, size(A, 2)).
    
    log : structure array
        Holds fields with logged information on routine - fields depend on
        subroutines used.
        
    ----------
    This algorithm computes Q and then sets B = Q.T @ A. Conceptually, we
    compute Q by using Algorithm 4.3 (see also Algorithm 4.4) from
        Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp.
        "Finding structure with randomness: Probabilistic algorithms for
        constructing approximate matrix decompositions."
        SIAM review 53.2 (2011): 217-288.
        (available at `arXiv <http://arxiv.org/abs/0909.4061>`_).
    The precise subspace iteration technique is similar to that of Algorithm
    3.3 from
         Bolong Zhang and Michael Mascagni.
         "Pass-Efficient Randomized LU Algorithms for Computing Low-Rank
         Matrix Approximation"
         arXiv:2002.07138 (2020).
    The main difference between this implementation and Zhang and
    Mascagni's Algorithm 3.3: we use QR decompositions where they use LU
    decompositions. 
    Additionally, using an alternative sketching scheme from
    '../rangefinders' would allow 0 steps or 1 step of
    subspace iteration, where Zhang and Mascagni's Algorithm
    implementation requires >= 2 steps.
    
    ----------
    Important note:
    before calling this routine, use:
    addpath('../utils') - for MatrlaRandStream.m
    addpath('../randgefinders') - for different versions of rangefinders.
%}

    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        log.status = 'Optional parameter for logging detailed information has not been passed.'; 
    else
        log_present = 1;
        logging.depth = logging.depth - 1;
    end

    s = MarlaRandStream(s);
    % Sketch construction stage - alternative options are available in 
    %'../rangefinders'.
    [Q, log] = rf1(A, k, p, s, logging);
    % Stage of computing B deterministically. 
    B = Q' * A;

    % Relative error computation
    if logging == 2 && log_present
        log.A_fro_norm = norm(A, 'fro');
        log.absolute_error = norm(A - Q * B, 'fro');
    end
end
