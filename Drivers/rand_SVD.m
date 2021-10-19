function [U, S, V] = rand_SVD(A, k, s, p)

    [Q, B] = rand_QB(A, k, s, p);
    [U_, S, V] = svd(B);
    U = Q * U_;
end