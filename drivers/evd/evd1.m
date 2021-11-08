% This function isn't fully implemented.
function[V, lambda_matrix] = evd1(A, k, s, p, tol, block_size)
%{

    Rely on a rangefinder to obtain the matrix Q for the decomposition
    A \approx Q B. Once we have Q, we construct B = Q' * A and return
    (Q, B). This function is agnostic to the implementation of the
    rangefinder: it might build a rank-k matrix Q all at once or construct
    successively larger matrices Q by an iterative process. We make no
    assumptions on the rangefinder's termination criteria beyond those
    listed below.
    Parameters
    ----------
    A : matrix
        Data matrix to approximate. A must be an n × n Hermitian matrix.
    k : int
        Target rank for the approximation of A: 0 < k < min(size(A)).
        This parameter includes any oversampling. For example, if you
        want to be near the optimal (Eckhart-Young) error for a rank 20
        approximation of A, then you might want to set k=25.
    s  : int
        Auxiliary parameter for the QBDecomposer or RangeFinder.
    p : int
        Number of power iterations used. Using this algorithm version
        implies increase of the number of passes over matrix A by 2 with each
        power iteration. 
    tol : float
        Target accuracy for the oversampled approximation of A: 0 < tol < np.inf.
        This parameter inherits from the QBDecomposer or RangeFinder class.

    Returns
    -------
    V : matrix
        Has size (size(A, 1), k). Columns are orthonormal.
    lambda_matrix : matrix
        Has size (k,k). lambda_matrix is a diagonal matrix storing the eigenvalues of 
        the matrix A.
%}
        % Using a version of QB algorithm. Alternative versions may be found in
        % '../Comps/QB'. 
        addpath('../Comps/QB');
        [Q, B] = rand_qb(A, k + s, p);
        % B = Q' * A is necessary.
        C = B * Q;
        % d = number of columns in Q, d ≤ k + s
        d = size(Q, 2);
        if d > (k + s)
            disp('This implementation has the dimension of Q matrix <= k + s.');
        end
        [U, lambda_matrix] = eig(C, 'vector');
        r = min(k, d);
        I = ismember(lambda_matrix, maxk(lambda_matrix, r));
        % Indices of r largest components of |λ|.
        U = U(:, I);
        lambda_matrix = lambda_matrix(I, :); 
        V = Q * U;
        lambda_matrix = diag(lambda_matrix);
end