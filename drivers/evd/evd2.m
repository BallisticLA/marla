function[V, lambda_matrix] = evd2(A, k, s, p)
%{
    Once we have Omega, we construct Q = A * Omega and return Cholesky decomposition of 
    Omega' * Q. This function is agnostic to the implementation of the
    sketching: it might has different ways of row sketching. We make no
    assumptions on the RowSketcher's termination criteria beyond those
    listed below.
    Return the eigen decomposition matrices (V, lambda_matrix),
    based on a row sketched version of A.
    Use a Gaussian sketching matrix and pass over A a total of
    num_passes times.
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

    Returns
    -------
    V : matrix
        Has size (size(A, 1), k). Columns are orthonormal.
    lambda_matrix : matrix
        Has size (k,k). lambda_matrix is a diagonal matrix storing the eigenvalues of 
        the matrix A.
%}
        % Sketch construction stage - alternative options are available in 
        %'../rangefinders'.
        class_A = class(A);
        [m, n] = size(A);
        % By default, a Gaussian random sketching matrix is used.
        % Alternative choices are present in '../Sketching_Operators'
        Omega = randn(n, k + s, class_A);
        [Q, ~] = qr(A * Omega, 0);
    
        for j = 1 : p
            [Q, ~] = qr(A' * Q, 0);
            [Q, ~] = qr(A * Q, 0);
        end

        nu = sqrt(m) * eps('double') * norm(Q, 'fro');
        % A temporary regularization parameter.
        Q = Q + nu * Omega;
        R = chol(Omega' * Q, 'lower');
        % R is upper-triangular and R' * R = Omega' * Q = Omega' * (A + nu
        % * I) * Omega
        B = Q  / R';
        B = (R \ Q')';
        % B has n rows and k + s columns.
        [V, Sigma_matrix, ~] = svd(B);
        
        comp_list = [k];
        for i = 1:min(k, n)
            if Sigma_matrix(i+1) ^ 2 <= nu
                comp_list = [comp_list i]; %#ok<AGROW>
            end
        end
        % comp_list constracuts the union from which we drop components next.        
        r = min(comp_list); 
        % drop components that relied on regularization
        lambda_matrix = ((Sigma_matrix(1:r)) ^ 2) - nu;
        V = V(:, 1:r);
        lambda_matrix = diag(lambda_matrix);
end