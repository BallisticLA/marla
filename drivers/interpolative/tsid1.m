function[Z, Is, X, Js, log] = tsid1(A, k, over, p, seed, logging)
%{
    Computes double-sided Interpolative Decomposition of matrix A.
    Rank-k approximation of A is then present as:
    A ~= Z * (A(Is(1 : k), Js(1 : k))) * X.

    Obtain a one-sided ID by any means, then deterministically extend
    to a two-sided ID.
    Using OSID1 would make this a "Sketch + QRCP" approach to double ID,
    as described in Voronin & Martinsson, 2016, Sections 2.4 and 4.
    
    ----------
    A : matrix
        Data matrix to approximate. Said to be of size m by n.

    k : int

        The returned approximation will be truncated to rank k.
    over : int
        Oversampling parameter.

    p : int
        Number of power iterations used. Using this algorithm version
        implies increase of the number of passes over matrix A by 2 with each
        power iteration.

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
    ----------
    Z : matrix
        Size m by k.

    Is : vector
        Size k by 1.
        Values represent select rows of matrix A.

    X : matrix
        Size k by n.

    Js : vector
        Size k by 1.
        Values represent select columns of matrix A.

    log : structure array
         Holds fields with logged information on routine - fields depend on
         subroutines used.

    Important note: 
    ---------------
    Before calling this routine, use:   
    addpath('../../utils') - for MatrlaRandStream.m
    addpath('../../comps/interpolative/');
%}
    seed = MarlaRandStream(seed);
    % Relies on a computational routine qrcp_osid
    if size(A, 1) > size(A, 2)
        [X, Js, log] = osid1(A, k, over, p, 1, seed, logging);
        [Z, Is] = qrcp_osid(A(:, Js), k, 0);
    else
        [Z, Is, log] = osid1(A, k, over, p, 0, seed, logging);
        [X, Js] = qrcp_osid(A(Is, :), k, 1);
    end
end