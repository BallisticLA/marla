function [Q, B] = rand_QB_SP(A, k)
%{
"Single-pass" version of QB algorithm.
Return matrices (Q, B) from a rank-k QB factorization of A.
----------
A : matrix
    Data matrix to approximate.
k : int
    Target rank for the approximation of A: 0 < k < min(size(A)).
    This parameter includes any oversampling. For example, if you
    want to be near the optimal (Eckhart-Young) error for a rank 20
    approximation of A, then you might want to set k=25.
Returns
-------
Q : matrix
    Has size (size(A, 1), k). Columns are orthonormal.
B : matrix
    Has shape (k, size(A, 2)).
----------
Avoids using power iterations by definition. 
Requires computing two sketches - in the provided implementation,
this stage is not done in a single pass, although it can be.
References
----------
Initial idea stated in section 5.5 of  Halko, Nathan, Per-Gunnar
Martinsson, and Joel A. Tropp's paper. Further developed in YGL:2018.
%}
    class_A = class(A);
    [m, n] = size(A);
    % Sketch construction stage - two sketches are constructed here. 
    Omega = randn(n, k, class_A);
    Omega_ = randn(m, k, class_A);
    Y = A * Omega;
    Y_ = A' * Omega_;
    [Q, ~] = qr(Y, 0);
    [Q_, ~] = qr(Y_, 0);
    % Stage of computing B deterministically.
    B = (Omega_' * Q) \ (Y_' * Q_) * Q_';
end