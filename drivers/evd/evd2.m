function[V, lambda_matrix] = evd2(A, k, over, num_passes)
%
%     Return a matrix V and column vector lambda that define a positive 
%     semidefinite matrix "A_approx" through its eigen-decomposition:
% 
%         A_approx = V * diag(lambda) * V'.
% 
%     The function assumes A is symmetric positive semidefinite.
%     The columns of V are approximations of the dominant eigenvectors of A.
%     The entries of lambda are the corresponding approximate eigenvalues.
% 
%     The approximation is "fixed rank." The matrix V has d = min(k, rank(A))
%     columns, and there is no direct control over the approximation error
%     ||A_approx - A||. Increasing "num_passes" and "over" should result in
%     better approximations.
% 
%     Parameters
%     ----------
%     A : matrix
%         Data matrix to approximate. A must be an n × n Hermitian matrix.
% 
%     k : int
%         Target rank for the approximation of A: 0 < k < min(size(A)).
%         This parameter includes any oversampling. For example, if you
%         want to be near the optimal (Eckhart-Young) error for a rank 20
%         approximation of A, then you might want to set k=25.
% 
%     over : int
%         Perform internal calculations with a sketch of rank (k + over).
%         This is usually a small constant, e.g., 5 to 25. In some situations
%         it's useful to set over = k.
% 
%     num_passes : int
%         Total number of passes the algorithm is allowed over A.
%         We require num_passes >= 2, and usually we have num_passes <= 10.
%         Increasing this parameter is one way to obtain better
%         approximations, especially at lower ranks.
% 
%     Returns
%     -------
%     V : matrix
%         Has size (n, d), where d = min(k, rank(A)).
%         Columns are orthonormal.
% 
%     lambda : vector
%         Has size (d, 1), where d = min(k, rank(A)).
%
%    Notes 
%    -----
%    Sketch construction stage - alternative options are available in 
%    '../rangefinders'.
%
    class_A = class(A);
    [m, n] = size(A);
    % By default, a Gaussian random sketching matrix is used.
    % Alternative choices are present in '../Sketching_Operators'
    Omega = randn(n, k + over, class_A);
    [Q, ~] = qr(A * Omega, 0);
    
    for j = 1 : num_passes
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