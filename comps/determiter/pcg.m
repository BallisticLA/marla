function [x, y, resid_vec] = pcg(A, b, c, delta, tol, iter_lim, M, x0)
% Return an approximate solution to the saddle point system with data
%   
%       [ I  ,      A    ] * [y] = [b]
%       [ A' , -delta * I]   [x]   [c].
%
%   The matrix A is m-by-n and tall. The vector b has length m and 
%   the vector c has length n. delta is a nonnegative scalar.
%
%   This function implements Algorithm B.3. from the appendix of
%   the monograph "Conjugate Gradients without the agonizing pain."
%
%       Our preconditioner M relates to the preconditioner in the monograph
%       by M * M' = (M_{mono})^{-1}
%
%       Our matrix A relates to the matrix in the monograph by
%       A' * A + delta*eye(n) = A_{mono}.
%
%       Our vector b relates to the vector b in the monograph by
%       b_{mono} = A'*b - c.
%
    b1 = A'*b - c;
    r = b1 - (A'*(A*x0) + delta*x0);
    d = M * (M'*r);
    reg = delta > 0;
    
    x = x0;
    delta1_old = r' * d;
    delta1_new = delta1_old;
    rel_sq_tol = delta1_old * (tol^2);
    
    resid_vec = ones(iter_lim);
    i = 1;
    while i <= iter_lim && delta1_new > rel_sq_tol
        resid_vec(i) = delta1_old;
        q = A'*(A*d);
        if reg
            q = q + delta*d;
        end
        alpha = delta1_new / (d'*q);
        x = x + alpha*d;
        if mod(i,25) == 0
           r = b1 - A'*(A*x);
           if reg
               r = r - delta*x;
           end
        else
           r = r - alpha*q;
        end
        s = M*(M'*r);
        delta1_old = delta1_new;
        delta1_new = r'*s;
        beta = delta1_new / delta1_old;
        d = s + beta*d;
        i = i + 1;
    end
    resid_vec = resid_vec(1:i);
    y = b - A*x;
end

