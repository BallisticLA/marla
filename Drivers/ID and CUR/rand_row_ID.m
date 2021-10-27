function [J, V] = rand_row_ID(A, k, s, p) 
    %{
    Computes row Interpolative Decomposition of matrix A.
    Rank-k approximation of A is then present as:
    A ~= V*(A(J(1 : (k + s)),:).
    
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
    
    Does not use sketch construction versions from 
    '../../Comps/Sketch_Construction' due to subtle difference at the power 
    iterations stage.
 
    Algorithm to be revised. 

    Reference
    ----------
    Section 3.5 of https://arxiv.org/pdf/1502.05366.pdf - RSVDPACK notes.
    %}
    class_A = class(A);
    [m, n] = size(A);
    l = k + s;
    % By default, a Gaussian random sketching matrix is used.
    % Alternative choices are present in '../../Comps/Sketching_Operators'
    Omega = randn(n, l, class_A);
    Y = A * Omega;
    %Power Iterations
    for j = 1 : p
        [Y, ~] = qr(transpose(A) * Y, 0);
        [Y, ~] = qr(A * Y, 0);
    end
    
    [~, R, J] = qr(Y', 'vector');
    
    I = eye(m, m);

    V = I(:, J) * [eye(l, l) R(:, 1 : l) \ R(:, l + 1 : m)]';
end