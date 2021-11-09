function[V, lambda] = evd1(A, k, tol, over, num_passes, block_size)
% 
%     Return a matrix V and column vector lambda that define a symmetric matrix
%     "A_approx" through its eigen-decomposition:
% 
%         A_approx = V * diag(lambda) * V'.
% 
%     The function assumes A is symmetric.
%     The columns of V are approximations of the dominant eigenvectors of A.
%     The entries of lambda are the corresponding approximate eigenvalues.
% 
%     Parameters
%     ----------
%     A : matrix
%         Data matrix to approximate. Must be n × n and symmetric matrix.
% 
%     k : int
%         Maximum rank for the approximation of A: 0 < k <= n.
%         If you set tol=0, then the returned approximation will have rank
%         min(k, rank(A)).
% 
%     tol : float
%         Target accuracy for the oversampled approximation of A: 0 < tol < np.inf.
%         If you set k = n and over=0 then the returned approximation should
%         satisfy ||A - A_approx||_F <= tol.
% 
%     over : int
%         Perform internal calculations with a sketch of rank (k + over).
%         This is usually a small constant, e.g., 5 to 25. In some situations
%         it's useful to set over = k. It's valid to set over = 0.
% 
%     num_passes : int
%         Total number of passes the algorithm is allowed over A.
%         We require num_passes >= 2, and usually we have num_passes <= 10.
%         Increasing this parameter is one way to obtain better
%         approximations, especially at lower ranks.
% 
%     block_size : int
%         The approximation is built incrementally, updating the rank by
%         block_size at each iteration (with safeguards so we never exceed
%         the rank of A).
% 
%     Returns
%     -------
%     V : matrix
%         Has size (n, d), where d <= min(k, rank(A)).
%         Columns are orthonormal.
% 
%     lambda : vector
%         Has size (d,1), where d <= min(k, rank(A)).
% 
%    Notes 
%    -----
%    Using a version of QB algorithm. Alternative versions may be found in
%    '../comps/qb'. 
%
    addpath('../comps/qb');
    [Q, B] = rand_qb_b(A, block_size, tol / 2, k + over, num_passes);
    % B = Q' * A is necessary.
    C = B * Q;
    % d = number of columns in Q, d ≤ k + over
    d = size(Q, 2);
    r = min(k, d);
    [U, lambda] = eig(C, 'vector');
    alambda = abs(lambda);
    r = min(r, nnz(alambda > eps('double')));
    [~,I] = sort(-alambda);
    I = I(1:r);
    % Indices of r largest components of |λ|.
    U = U(:, I);
    lambda = lambda(I); 
    V = Q * U;
end
