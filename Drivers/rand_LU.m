function [P1, P2, L, U] = rand_LU(A, k, s, p)

    [Q, B] = rand_QB(A, k, s, p);
    [L1, U1, P2] = lu(B');
    Y = Q * U1';
    [L2, U2, P1] = lu(Y);
    L = L2;
    U = U2 * L1';
end