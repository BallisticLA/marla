function [U, S, V] = svd1(A, k, s, p, tol, block_size)
%{
Return U, S, V where, for some integer ell <= k,
    U is size(A, 1)-by-ell,
    S is a vector of length ell,
    V is ell-by-size(A, 2),
    so that
    A \approx U * diag(S) * V
Parameters
----------
A : matrix
    Data matrix to approximate
k : int
    The returned SVD will be truncated to at most rank k:
    0 < k <= min(size(A)). Setting k=min(size(A)) and s=0
    ensures ||A - U * diag(S) * V|| / ||A|| <= tol on exit. However,
    setting k=min(size(A)) may trivially return the SVD of
    A in some implementations.
s : int
    The randomized part of the algorithm uses k+s as the target rank;
    we require over >= 0 and k+s <= min(size(A)).
    In a conformant implementation, that part of the algorithm will
    never return a factorization of rank greater than k+over.
    Setting over > 0 will likely result in truncating the SVD obtained
    from the randomized part of the algorithm. If you want to control
    the truncation step yourself, then you should set over=0 and
    increase the value of k. E.g., a function call with over=5 and
    k=20 can avoid truncation by setting k=25 and over=0.
p : int
    Number of power iterations used. Using this algorithm version
    implies increase of the number of passes over matrix A by 2 with each
    power iteration. 
tol : float
    The target error used by the randomized part of the algorithm.
    When over = 0, this parameter controls ||A - U @ diag(s) @ Vh||.
    The precise meaning of "tol" when over > 0 is implementation
block_size : int
    The block size for a blocked QB algorithm. Add this many columns
    to Q at each iteration (except possibly the final iteration).
    dependent.

This implementation of rand_SVD uses versions of QB algorithm as its main 
computational routine. Tolerance and block size parameters are not required 
for some QB implementations.
%}
    % Using a version of QB algorithm. Alternative versions may be found in
    % '../comps/qb'. 
    addpath('../comps/qb');
    [Q, B] = rand_qb_b(A, block_size, tol, k, p);
    % Using a built-in function for computing an SVD. 
    [U, S, V] = svd(B);

    % Removing singular values below machine precision. 
    cutoff = find(diag(S) < eps(class(S)), 1);
    % Removing "oversampled" data. 
    if s > 0
        cutoff = min(k, cutoff);
    end
    if(~isempty(cutoff))
        U = U(:, 1:cutoff);
        S = S(1:cutoff, 1:cutoff);
        V = V(:, 1:cutoff);
    end
    
    % Adjusting matrix U.
    U = Q * U;
end
