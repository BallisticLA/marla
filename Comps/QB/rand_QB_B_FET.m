function [Q, B] = rand_QB_B_FET(A, block_size, tol, k, p)
    norm_A = norm(A, 'fro');
    
    if norm_A == 0
        fprintf('The input matrix is empty.');
        return
    end
    
    class_A = class(A);
    [m, n] = size(A);
    norm_B = 0;
    Q = zeros(m, 0, class_A);
    B = zeros(0, n, class_A);
    
    for i = 1 : (k / block_size)
    
        Omega_i = randn(n, block_size, class_A);
        [Q_i, ~] = qr((A * Omega_i) - (Q * (B * Omega_i)), 0);
        
        for j = 1 : p
            [Q_i, ~] = qr(A' * Q_i - B' * (Q' * Q_i), 0);
            [Q_i, ~] = qr(A * Q_i - Q * (B * Q_i), 0);
        end
        
        [Q_i, ~] = qr(Q_i - (Q * (Q' * Q_i)), 0);
        B_i = Q_i' * A;
        
        norm_B = hypot(norm_B, norm(B_i, 'fro'));
        prev_error = approximation_error;
        approximation_error = sqrt(abs(norm_A - norm_B) * (norm_A + norm_B)) / norm_A;
        
        %Handling the round-off error accumulation
        if (i > 1) && (approximation_error > prev_error)
            break
        end
        
        Q = [Q, Q_i]; %#ok<AGROW>
        B = [B; B_i]; %#ok<AGROW>
        
        if approximation_error < tol       
            break;
        end
    end
    
    if i == (k / block_size)
        fprintf('Approximation error = %f. Fail to converge within the specified tolerance.\n\n', approximation_error / norm_A);
    end
end