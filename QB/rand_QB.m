function [Q, B] = rand_QB(A, k, s, p)

    class_A = class(A);
    [~, n] = size(A);
    Omega = randn(n, k + s, class_A);
    [Q, ~] = qr(A * Omega, 0);
    
    % Power Iterations.
    for j = 1 : p
        [Q, ~] = qr(A' * Q, 0);
        [Q, ~] = qr(A * Q, 0);
    end
    
    B = Q' * A;
end