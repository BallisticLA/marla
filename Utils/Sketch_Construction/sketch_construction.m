function [Q] = sketch_construction(A, Omega, p)

    [Q, ~] = qr(A * Omega, 0);

    for j = 1 : p
        [Q, ~] = qr(A' * Q, 0);
        [Q, ~] = qr(A * Q, 0);
    end
end