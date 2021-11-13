function [Q, B] = rand_qb(A, k, p, s)
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
Returns
-------
Q : matrix
    Has size (size(A, 1), k). Columns are orthonormal.
B : matrix
    Has shape (k, size(A, 2)).
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
%}
    s = MarlaRandStream(s);
    % Sketch construction stage - alternative options are available in 
    %'../rangefinders'.
    class_A = class(A);
    [~, n] = size(A);
    % By default, a Gaussian random sketching matrix is used.
    % Alternative choices are present in '../Sketching_Operators'
    Omega = randn(s, n, k, class_A);
    [Q, ~] = qr(A * Omega, 0);

    for j = 1 : p
        [Q, ~] = qr(A' * Q, 0);
        [Q, ~] = qr(A * Q, 0);
    end
    % Stage of computing B deterministically. 
    B = Q' * A;
end
