function [Q, R] = rand_QR(A, k, s, p)

    [Q, B] = rand_QB(A, k, s, p);
    [Q_, R] = qr(B);
    Q = Q * Q_;
end