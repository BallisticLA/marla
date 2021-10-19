function [J, V] = rand_row_ID(A, k, s, p) 
    
    class_A = class(A);
    [m, n] = size(A);
    l = k + s;
    
    Omega = randn(n, l, class_A);
    
    Y = A * Omega;
    
    
    for j = 1 : p
        [Y, ~] = qr(transpose(A) * Y, 0);
        [Y, ~] = qr(A * Y, 0);
    end
    
    [~, R, J] = qr(Y', 'vector');
    
    I = eye(m, m);

    V = I(:, J) * [eye(l, l) R(:, 1 : l) \ R(:, l + 1 : m)]';
end