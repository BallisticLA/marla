function[V, lambda, log] = evd2(A, k, over, num_passes, s, logging)
%{
     Return a matrix V and column vector lambda that define a positive 
     semidefinite matrix "A_approx" through its eigen-decomposition:
 
         A_approx = V * diag(lambda) * V'.
 
     The function assumes A is symmetric positive semidefinite.
     The columns of V are approximations of the dominant eigenvectors of A.
     The entries of lambda are the corresponding approximate eigenvalues.
 
     The approximation is "fixed rank." The matrix V has d = min(k, rank(A))
     columns, and there is no direct control over the approximation error
     ||A_approx - A||. Increasing "num_passes" and "over" should result in
     better approximations.
 
     Parameters
     ----------
     A : matrix
         Data matrix to approximate. A must be an n × n Hermitian matrix.
 
     k : int
         Target rank for the approximation of A: 0 < k < min(size(A)).
         This parameter includes any oversampling. For example, if you
         want to be near the optimal (Eckhart-Young) error for a rank 20
         approximation of A, then you might want to set k=25.
 
     over : int
         Perform internal calculations with a sketch of rank (k + over).
         This is usually a small constant, e.g., 5 to 25. In some situations
         it's useful to set over = k.
 
     num_passes : int
         Total number of passes the algorithm is allowed over A.
         We require num_passes >= 2, and usually we have num_passes <= 10.
         Increasing this parameter is one way to obtain better
         approximations, especially at lower ranks.

     logging : struct array 
         Parameter for logging different levels of detailed information.
         Contains two fields:
         Depth - describes on how many subroutines levels should the info
         be logged (0 for none, 1 for just the main routine, etc).
         Span - describes which types of info should be logged (0 for none,
         1 for timings, 2 for timings + errors estimates, 3 for unusual 
         behaviors).
 
     Returns
     -------
     V : matrix
         Has size (n, d), where d = min(k, rank(A)).
         Columns are orthonormal.
 
     lambda : vector
         Has size (d, 1), where d = min(k, rank(A)).

     log : structure array
         Holds fields with logged information on routine - fields depend on
         subroutines used.

    Notes 
    -----
    Sketch construction stage - alternative options are available in 
    '../rangefinders'.
    Before calling this routine, use:
    addpath('../utils') - for MatrlaRandStream.m
%}
    if logging.depth == 0 || logging.span == 0
        log_present = 0;
        %disp('Optional parameter for logging detailed information has not been passed.'); 
    else
        log_present = 1;
        logging.depth = logging.depth - 1;
    end

    s = MarlaRandStream(s);
    class_A = class(A);
    n = size(A, 1);
    assert(k < n);
    % By default, a Gaussian random sketching matrix is used.
    % Alternative choices are present in '../Sketching_Operators'
    if log_present, tic, end
    Omega = randn(s, n, k + over, class_A);
    [Q, ~] = qr(A * Omega, 0);
    if log_present, log.t_sketch = toc; end
    % Power Iterations
    if log_present, tic, end
    for j = 1 : num_passes
        [Q, ~] = qr(A' * Q, 0);
        [Q, ~] = qr(A * Q, 0);
    end
    if log_present, log.t_power_iter = toc; end
    % Main part of the routine
    if log_present, tic, end
    Y = A * Q;
    nu = sqrt(n) * eps('double') * norm(Y, 'fro');
    % A temporary regularization parameter.
    Y = Y + nu * Q;
    R = chol(Q' * Y, 'lower');
    B = Y  / R';
    % B has n rows and k + s columns.
    [V, Sigma_matrix, ~] = svd(B, 'econ');
    sigma = diag(Sigma_matrix);
    
    comp_list = [k];
    for i = 1:(k-1)
        if sigma(i+1) ^ 2 <= nu
            comp_list = [comp_list i]; %#ok<AGROW>
        end
    end
    % comp_list constracuts the union from which we drop components next.        
    r = min(comp_list); 
    % drop components that relied on regularization
    lambda = ((sigma(1:r)).^ 2) - nu;
    V = V(:, 1:r);
    if logging, log.t_main = toc; end

    % Relative error computation
    if logging.span >= 2 && log_present
        log.A_fro_norm = norm(A, 'fro');
        log.absolute_error = norm(A - V * diag(lambda) * V', 'fro');
    end
end