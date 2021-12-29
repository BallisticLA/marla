function[Js, U, Is, log] = cur1(A, k, over, p, seed, logging)
%{
    Computes CUR Decomposition of matrix A.
    Relies on row ID algorithm.

    Rank-k approximation of A is then present as:
    A ~= Js * U * Is.

    Input
    ----------
    A : matrix
        Data matrix to approximate

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
    Js : matrix
        
    U : matrix
    
    Is : matrix

    log : structure array
         Holds fields with logged information on routine - fields depend on
         subroutines used.

    Reference
    ----------
    Section 2.5 of https://arxiv.org/pdf/1502.05366.pdf - RSVDPACK notes.
%}  
    seed = MarlaRandStream(seed);
    if size(A, 1) > size(A, 2)
        [X, Js, log] = osid1(A, k, over, p, 1, seed, logging);
        [~, ~, Is] = qr(A(:, Js)', 0);
        Is = Is(:);  % to column vector
        Is = Is(1 : k);
        U = X / A(Is, :);
    else
        [Z, Is, log] = osid1(A, k, over, p, 0, seed, logging);
        [~, ~, Js] = qr(A(Is, :), 0);
        Js = Js(:);  % to row vector
        Js = Js(1:k);
        U = A(:, Js) \ Z;
    end
end