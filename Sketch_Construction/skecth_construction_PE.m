function [V] = skecth_construction_PE(A, k, s, power)

    class_A = class(A);
    [m, n] = size(A);
    l = k + s;
    v = 2 * power + 1;

    % Odd number of passes over A.
    if(mod(v, 2) == 0)
        
        Omega = randn(m, l, class_A);
        
        if (v > 2)
            [V, ~] = lu(A' * Omega);
        else
            [V, ~] = qr(A' * Omega);
        end
    % Even number of passes over A.
    else
        V = randn(n, l, class_A);
    end

    for i = 1 : power

        [V, ~] = lu(A * V);

        if i == power
            [V, ~] = qr(A' * V);
        else
            [V, ~] = lu(A' * V);
        end
    end
end