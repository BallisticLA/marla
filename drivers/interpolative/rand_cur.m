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
    [J_r, ~] = rand_row_ID(A, k, s, p);
    R = A(J_r(1 : (k + s)), :);
    
    [J_c, V_c] = rand_row_ID(R', k, s, p);
    C = A(:, J_c(:, 1 : (k + s)));
    
    U = V_c' / R;
end