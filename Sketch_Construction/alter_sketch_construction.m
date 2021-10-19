function [Q] = alter_sketch_construction(A, Omega, p)

    [Q, ~] = qr(A * Omega, 0);

    for j = 1 : p
    
        [Q, ~] = lu(A' * Q);
        
        if j < power 
            [Q, ~] = lu(A * Q);
        else
            [Q, ~] = qr(A * Q, 0);
        end
    end
end