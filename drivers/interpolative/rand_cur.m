function [C, U, R] = rand_cur(A, k, s, p) 
    %{
    Computes CUR Decomposition of matrix A.
    Relies on row ID algorithm.

    Rank-k approximation of A is then present as:
    A ~= C * U * R.

    Parameters
    ----------
    A : matrix
        Data matrix to approximate
    k : int
        The returned approximation will be truncated to rank k.
    s : int
        Oversampling parameter.
    p : int
        Number of power iterations used. Using this algorithm version
        implies increase of the number of passes over matrix A by 2 with each
        power iteration.

    Reference
    ----------
    Section 2.5 of https://arxiv.org/pdf/1502.05366.pdf - RSVDPACK notes.
    %}
    addpath('../../comps/interpolative/');
    Js = rocs1(A, k, s, p, 0);
    R = A(Js(1 : k), :);

    % By default, uses osid1. 
    [Z, Is] = osid1(R', k, s, p, 0);
    C = A(:, Is(:, 1 : k));
    
    U = Z' / R;
end