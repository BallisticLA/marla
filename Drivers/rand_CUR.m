function [C, U, R] = rand_CUR(A, k, s, p) 
    
    [J_r, ~] = rand_row_ID(A, k, s, p);
    R = A(J_r(1 : (k + s)), :);
    
    [J_c, V_c] = rand_row_ID(R', k, s, p);
    C = A(:, J_c(:, 1 : (k + s)));
    
    U = V_c' / R;
    
end