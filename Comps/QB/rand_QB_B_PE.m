function [Q, B] = rand_QB_B_PE(A, block_size, tol, k, p)
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

We start by obtaining a sketching matrix of shape (size(A, 1), k),
using p steps of subspace iteration on a random Gaussian
matrix with k columns. Then we perform two more passes over A before
beginning iterative construction of (Q, B). Each iteration adds at most
"block_size" columns to Q and rows to B.
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
References
----------
This implements a variant of [YGL:2018, Algorithm 4].
%}
    norm_A = norm(A, 'fro');
    % Early termination check on an empty input. 
    if norm_A == 0
        fprintf('The input matrix is empty.');
        return
    end
    class_A = class(A);
    [m, n] = size(A);
    norm_B = 0;
    % Sketch construction stage.
    Omega = randn(n, k, class_A);
    for  j = 1 : p
        [G, ~] = qr(A * Omega, 0);
        [Omega, ~] = qr(A' * G, 0);
    end
    G = A * Omega;
    H = A' * G;
    % Pre-initialization of output matrices. 
    Q = zeros(m, 0, class_A);
    B = zeros(0, n, class_A);
    % Iterative stage.
    curr_idx = 1;
    while curr_idx < k
        %Avoiding exceeding array bounds 
        if(block_size + curr_idx - 1 > k)
            block_size = l - curr_idx + 1;
        end
    
        Omega_i = Omega(:, curr_idx : curr_idx + block_size - 1);
        
        Temp = B * Omega_i;
        
        Y_i = G(:, curr_idx : curr_idx + block_size - 1) - (Q * Temp);
        [Q_i, R_] = qr(Y_i, 0);
        
        [Q_i, R_i] = qr(Q_i - (Q * (Q' * Q_i)), 0);
        R_i = R_i * R_;
        % Deterministic computation of B_i.
        B_i = transpose(R_i) \ ((H(:, curr_idx : curr_idx + block_size - 1))' - (Y_i' * Q) * B - (Temp' * B));
        % Approximation error check.
        norm_B = hypot(norm_B, norm(B_i, 'fro'));
        prev_error = approximation_error;
        approximation_error = sqrt(abs(norm_A - norm_B) * (norm_A + norm_B)) / norm_A;
        %Handling the round-off error accumulation
        if (curr_idx > 1) && (approximation_error > prev_error)
            break
        end
        % Output update. 
        Q = [Q, Q_i]; %#ok<AGROW>
        B = [B; B_i]; %#ok<AGROW>
        % Condition of reaching tolerance. 
        if approximation_error < tol       
            break;
        end 
        curr_idx = curr_idx + block_size;
    end
end