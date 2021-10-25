function [Q, B] = rand_QB_B_PE_FET(A, block_size, tol, k, p)
    
    norm_A = norm(A, 'fro');
    if norm_A == 0
        fprintf('The input matrix is empty.');
        return
    end

    class_A = class(A);
    [m, n] = size(A);
    norm_B = 0;
    
    Omega = randn(n, k, class_A);
    
    Q = zeros(m, 0, class_A);
    B = zeros(0, n, class_A);
    
    for  j = 1 : p
        [G, ~] = qr(A * Omega, 0);
        [Omega, ~] = qr(A' * G, 0);
    end
    
    G = A * Omega;
    H = A' * G;
    
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

        B_i = transpose(R_i) \ ((H(:, curr_idx : curr_idx + block_size - 1))' - (Y_i' * Q) * B - (Temp' * B));
        
        norm_B = hypot(norm_B, norm(B_i, 'fro'));
        prev_error = approximation_error;
        approximation_error = sqrt(abs(norm_A - norm_B) * (norm_A + norm_B)) / norm_A;
        
        %Handling the round-off error accumulation
        if (curr_idx > 1) && (approximation_error > prev_error)
            break
        end
        
        Q = [Q, Q_i]; %#ok<AGROW>
        B = [B; B_i]; %#ok<AGROW>
        
        if approximation_error < tol       
            break;
        end
        
        curr_idx = curr_idx + block_size;
    end
    
    if curr_idx == k
        fprintf('Approximation error = %f. Fail to converge within the specified tolerance.\n\n', approximation_error / norm_A);
    end
end