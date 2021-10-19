function [Q, B] = rand_QB_B_FR(A, k, s, p, block_size)
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

    Q = zeros(m, 0, class_A);
    B = zeros(0, n, class_A);
    
    for i = 1 : (l / block_size)
    
        Omega_i = randn(n, block_size, class_A);
        [Q_i, ~] = qr((A * Omega_i) - (Q * (B * Omega_i)), 0);
        
        for j = 1 : p
            [Q_i, ~] = qr(A' * Q_i - B' * (Q' * Q_i), 0);
            [Q_i, ~] = qr(A * Q_i - Q * (B * Q_i), 0);
        end
        
        [Q_i, ~] = qr(Q_i - (Q * (Q' * Q_i)), 0);
        B_i = Q_i' * A;
        
        Q = [Q, Q_i]; %#ok<AGROW>
        B = [B; B_i]; %#ok<AGROW>
        
        %In this case, keeping track of an approximation error is only used
        %for controlling the rouns-off error accumulation
        norm_B = hypot(norm_B, norm(B_i, 'fro'));
        prev_error = approximation_error;
        approximation_error = sqrt(abs(norm_A - norm_B) * (norm_A+norm_B)) / norm_A;
        
        %Handling the round-off error accumulation
        if (i > 1) && (approximation_error > prev_error)
            Q(:, end - block_size + 1 : end) = [];
            B(end - block_size + 1 : end, :) = [];
            break
        end
    end
end