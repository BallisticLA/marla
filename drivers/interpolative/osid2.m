function[Out1, Out2] = osid2(A, k, s, p, axis)
%{
    Sketch + (QRCP skeleton) + (least squares) approach to ID
    See Dong & Martinsson, 2021.

    Column ID returns an approximation of the form:
    A ~= (A(:, Out1(1 : k))) * Out2

    Row ID returns an approximation of the form:
    A ~= Out2 * (A(Out1(1 : k), :))
%}
    % Uses rf1 as a default rangefinder, but alternatives are available. 
    addpath('../../comps/rangefinders/');
    if axis == 0
        % Row ID
        S = rs1(A, k + s, p);
        Y = A * S;
        [~, ~, I] = qr(Y', 'vector');
        Out2 = I(1 : k);
        Out2 = Out2(:);
        Out1 = A / A(Out2, :);
        %return X, Is
    elseif axis == 1
        % Column ID
        S = rs1(A', k + s, p)';
        Y = S * A;
        [~, ~, J] = qr(Y, 'vector');
        Out2 = J(1 : k);
        Out2 = Out2(:);
        Out1 = A(:, Out2) \ A; 
        % Return Z, Js
    end
end