function[Out1, Out2, log] = osid1(A, k, over, p, axis, seed, logging)
%{
    Sketch + QRCP approach to ID
    See Voronin & Martinsson, 2016, Section 5.1.

    Column ID returns an approximation of the form:
    A ~= (A(:, Out1(1 : k))) * Out2

    Row ID returns an approximation of the form:
    A ~= Out2 * (A(Out1(1 : k), :))

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

    axis : int
        Taking values 0 or 1, choose row or column ID, respectively.

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
    Out1 : matrix
        
    Out2 : matrix

     log : structure array
         Holds fields with logged information on routine - fields depend on
         subroutines used.

    Important note: 
    ---------------
    Before calling this routine, use:
    addpath('../../comps/interpolative/'); 
    addpath('../../comps/rangefinders/');
%}
    seed = MarlaRandStream(seed);
    % Relies on a computational routine qrcp_osid
    if axis == 0
        % Row ID
        [S, log] = rs1(A, k + over, p, seed, logging);
        Y = A * S;
        [Out1, Out2] = qrcp_osid(Y, k, 0);
    elseif axis == 1
        % Column ID
        [S, log] = rs1(A', k + over, p, seed, logging);
        Y = S' * A;
        [Out1, Out2] = qrcp_osid(Y, k, 1);
    end
end