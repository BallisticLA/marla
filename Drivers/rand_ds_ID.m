function [J_r, J_c, V_r, V_c] = rand_ds_ID(A, k, s, p) 
    
    [J_r, V_r] = rand_row_ID(A, k, s, p);
    [J_c, V_c] = rand_row_ID(A(J_r(1 : (k + s)), :)', k, s, p);
    
end