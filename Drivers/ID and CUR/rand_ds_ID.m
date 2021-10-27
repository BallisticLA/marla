function [J_r, J_c, V_r, V_c] = rand_ds_ID(A, k, s, p) 
    %{
    Computes double-sided Interpolative Decomposition of matrix A.
    Rank-k approximation of A is then present as:
    A ~= V_r * (A(J_r(1 : (k + s)), J_c(1 : (k + s)))) * V_c'.

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
    % Row ID.
    [J_r, V_r] = rand_row_ID(A, k, s, p);
    % Column ID.
    [J_c, V_c] = rand_row_ID(A(J_r(1 : (k + s)), :)', k, s, p);  
end