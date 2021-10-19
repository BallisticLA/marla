function [Q, B] = rand_QB_SP(A, k, s)

    class_A = class(A);
    [m, n] = size(A);
    
    Omega = randn(n, k + s, class_A);
    Omega_ = randn(m, k + s, class_A);
    
    Y = A * Omega;
    Y_ = A' * Omega_;
    
    [Q, ~] = qr(Y, 0);
    [Q_, ~] = qr(Y_, 0);   
    B = (Omega_' * Q) \ (Y_' * Q_) * Q_';
end