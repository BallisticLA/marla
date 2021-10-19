function [Q, B] = rand_QB_B_PE_FR(A, k, s, p, block_size) 
    norm_A = norm(A, 'fro');
    
    if norm_A == 0
        fprintf('The input matrix is empty.');
        return
    end

    class_A = class(A);
    [m, n] = size(A);
    
    l = k + s;
    approximation_error = 0;
    norm_B = 0;
    
    Omega = randn(n, l, class_A);
    
    Q = zeros(m, 0, class_A);
    B = zeros(0, n, class_A);
    
    for  j = 1 : p
        [G, ~] = qr(A * Omega, 0);
        [Omega, ~] = qr(A' * G, 0);
    end
    
    G = A * Omega;
    H = A' * G;
    
    curr_idx = 1;
    while curr_idx < l
        
        %Avoiding exceeding array bounds 
        if(block_size + curr_idx - 1 > l)
            block_size = l - curr_idx + 1;
        end
        
        Omega_i = Omega(:, curr_idx : curr_idx + block_size - 1);
        
        Temp = B * Omega_i;
        
        Y_i = G(:, curr_idx : curr_idx + block_size - 1) - (Q * Temp);
        [Q_i, R_] = qr(Y_i, 0);
        
        [Q_i, R_i] = qr(Q_i - (Q * (Q' * Q_i)), 0);
        R_i = R_i * R_;

        B_i = transpose(R_i) \ ((H(:, curr_idx : curr_idx + block_size - 1))' - (Y_i' * Q) * B - (Temp' * B));

        Q = [Q, Q_i]; %#ok<AGROW>
        B = [B; B_i]; %#ok<AGROW>
        
        %In this case, keeping track of an approximation error is only used
        %for controlling the rouns-off error accumulation
        norm_B = hypot(norm_B, norm(B_i, 'fro'));
        prev_error = approximation_error;
        approximation_error = sqrt(abs(norm_A - norm_B) * (norm_A+norm_B)) / norm_A;
        
        %Handling the round-off error accumulation
        if (curr_idx > 1) && (approximation_error > prev_error)
            Q(:, end - block_size + 1 : end) = [];
            B(end - block_size + 1 : end, :) = [];
            break
        end
        
        curr_idx = curr_idx + block_size;
        
    end
end