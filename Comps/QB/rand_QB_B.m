function [Q, B] = rand_qb_b(A, block_size, tol, k, p)
%{
Iteratively build an approximate QB factorization of A,
which terminates once one of the following conditions
is satisfied
    (1)  || A - Q B ||_Fro / ||A|| <= tol
or
    (2) Q has k columns.
or
    (2) || A - Q B ||_Fro / ||A|| has increased in comparison to previous
    iteration's result (unlikely case). 
Each iteration involves sketching A from the right by a sketching
matrix with "block_size" columns, and using "p" power iterations. 
Parameters
----------
A : matrix
    Data matrix to approximate.
block_size : int
    The block size in this blocked QB algorithm. Add this many columns
    to Q at each iteration (except possibly the final iteration).
tol : float
    Terminate if ||A - Q B||_Fro / ||A|| <= tol. Setting k = min(size(A)) is a
    valid way of ensuring ||A - Q B||_Fro / ||A|| <= tol on exit.
k : int
    Terminate if size(Q, 2) == k. Assuming k < rank(A), setting tol=0 is a
    valid way of ensuring size(Q, 2) == k on exit.
p : int
    Number of power iterations used. Using this algorithm version
    implies increase of the number of passes over matrix A by 2 with each
    power iteration. 
Returns
-------
Q : matrix
    Has the same number of rows of A, and orthonormal columns.
B : matrix
    Has the same number of columns of A.
Notes
-----
The number of columns in Q increase by "block_size" at each iteration, unless
that would bring size(Q, 2) > k. In that case, the final iteration only
adds enough columns to Q so that size(Q, 2) == k.
We perform p steps of subspace iteration for each
block of the QB factorization. We stabilize subspace iteration with
QR factorization at each step.
References
----------
Algorithm 2 from YGL:2018.
%}
    norm_A = norm(A, 'fro');
    % Early termination check on an empty input. 
    if norm_A == 0
        fprintf('The input matrix is empty.');
        return
    end
    % Setting initial error to zero.
    approximation_error = 0;
    class_A = class(A);
    [m, n] = size(A);
    norm_B = 0;
    % Pre-initialization of output matrices. 
    Q = zeros(m, 0, class_A);
    B = zeros(0, n, class_A);
    % Iterative stage.
    for i = 1 : (k / block_size)
        % Consstructiong a sketch for current iteration. 
        Omega_i = randn(n, block_size, class_A);
        [Q_i, ~] = qr((A * Omega_i) - (Q * (B * Omega_i)), 0);
        % Power iterations. 
        for j = 1 : p
            [Q_i, ~] = qr(A' * Q_i - B' * (Q' * Q_i), 0);
            [Q_i, ~] = qr(A * Q_i - Q * (B * Q_i), 0);
        end
        % Ensuring orthogonalization of Q. 
        [Q_i, ~] = qr(Q_i - (Q * (Q' * Q_i)), 0);
        % Deterministic computation of B_i.
        B_i = Q_i' * A;
        % Approximation error check.
        norm_B = hypot(norm_B, norm(B_i, 'fro'));
        prev_error = approximation_error;
        approximation_error = sqrt(abs(norm_A - norm_B) * (norm_A + norm_B)) / norm_A; 
        % Handling the round-off error accumulation.
        if (i > 1) && (approximation_error > prev_error)
            break
        end
        % Output update. 
        Q = [Q, Q_i]; %#ok<AGROW>
        B = [B; B_i]; %#ok<AGROW>
        % Condition of reaching tolerance. 
        if approximation_error < tol       
            break;
        end
    end
end
