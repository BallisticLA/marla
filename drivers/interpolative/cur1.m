function[Js, U, Is] = cur1(A, k, s, p)
%{
    Computes CUR Decomposition of matrix A.
    Relies on row ID algorithm.

    Rank-k approximation of A is then present as:
    A ~= Js * U * Is.

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

    if size(A, 1) > size(A, 2)
        [X, Js] = osid1(A, k, s, p, 1);
        [~, ~, Is] = qr(A(:, Js)', 0);
        Is = Is(:);  % to column vector
        Is = Is(1 : k);
        U = X / A(Is, :);
    else
        [Z, Is] = osid1(A, k, s, p, 0);
        [~, ~, Js] = qr(A(Is, :), 0);
        Js = Js(:);  % to row vector
        Js = Js(1:k);
        U = A(:, Js) \ Z;
    end
end