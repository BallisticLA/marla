function [U, S, V, log] = svd1(A, k, tol, over, p, block_size, seed, logging)
%{
    Return U, S, V where, for some integer ell <= k,
        U is size(A, 1)-by-ell,
        S is a vector of length ell,
        V is ell-by-size(A, 2),
        so that
        A \approx U * diag(S) * V

    Input
    -----
    A : matrix
        Data matrix to approximate

    k : int
        The returned SVD will be truncated to at most rank k:
        0 < k <= min(size(A)). Setting k=min(size(A)) and s=0
        ensures ||A - U * diag(S) * V|| / ||A|| <= tol on exit. However,
        setting k=min(size(A)) may trivially return the SVD of
        A in some implementations.

    over : int
        The randomized part of the algorithm uses k+s as the target rank;
        we require over >= 0 and k+over <= min(size(A)).
        In a conformant implementation, that part of the algorithm will
        never return a factorization of rank greater than k+over.
        Setting over > 0 will likely result in truncating the SVD obtained
        from the randomized part of the algorithm. If you want to control
        the truncation step yourself, then you should set over=0 and
        increase the value of k. E.g., a function call with over=5 and
        k=20 can avoid truncation by setting k=25 and over=0.

    p : int
        Number of power iterations used. Using this algorithm version
        implies increase of the number of passes over matrix A by 2 with each
        power iteration. 

    tol : float
        The target error used by the randomized part of the algorithm.
        When over = 0, this parameter controls ||A - U @ diag(s) @ Vh||.
        The precise meaning of "tol" when over > 0 is implementation.

    block_size : int
        The block size for a blocked QB algorithm. Add this many columns
        to Q at each iteration (except possibly the final iteration).
        dependent.

    seed: int or RandStream
         Seed for rand stream.

    logging : struct array 
        Parameter for logging different levels of detailed information.
        Contains two fields:
        Depth - describes on how many subroutines levels should the info
        be logged (0 for none, 1 for just the main routine, etc).
        Span - describes which types of info should be logged (0 for none,
        1 for timings, 2 for timings + errors estimates, 3 for unusual 
        behaviors).

     Output
     ------
     U : matrix
         Orthonormal
 
     S : matrix
         Diagonal matrix of singular values

     V : matrix
         Orthonormal

     log : structure array
         Holds fields with logged information on routine - fields depend on
         subroutines used.
    
    This implementation of rand_SVD uses versions of QB algorithm as its main 
    computational routine. Tolerance and block size parameters are not required 
    for some QB implementations.
    
    Important note:
    ---------------
    before calling this routine, use:
    addpath('../utils') - for MatrlaRandStream.m
    addpath('../comps/qb') - for different versions of QB algorithm.
%}
    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        log.status = 'Optional parameter for logging detailed information has not been passed.'; 
    else
        log_present = 1;
        logging.depth = logging.depth - 1;
    end

    seed = MarlaRandStream(seed);
    % Using a version of QB algorithm. Alternative versions may be found in
    % '../comps/qb'. 
    if log_present, tic, end
    [Q, B, log] = rand_qb_b(A, block_size, tol, k, p, seed, logging);
    if log_present, log.t_qb = toc; end
    % Using a built-in function for computing an SVD. 
    if log_present, tic, end
    [U, S, V] = svd(B, 'econ');
    if log_present, log.svd = toc; end

    % Removing singular values below machine precision. 
    cutoff = find(diag(S) < eps(class(S)), 1);
    % Removing "oversampled" data. 
    if over > 0
        cutoff = min(k, cutoff);
    end
    if(~isempty(cutoff))
        U = U(:, 1:cutoff);
        S = S(1:cutoff, 1:cutoff);
        V = V(:, 1:cutoff);
    end
    
    % Adjusting matrix U.
    U = Q * U;
end
