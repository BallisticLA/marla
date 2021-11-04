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
        Y = rf1(A, k + s, p);
        [~, ~, I] = qr(Y', 'vector');
        Out1 = I(1 : k);
        Out2 = A / A(Out1, :);
        %return X, Is
    elseif axis == 1
        % Column ID
        Y = rf1(A', k + s, p);
        [~, ~, J] = qr(Y, 'vector');
        Out1 = J(1 : k);
        Out2 = A(:, Out1) \ A; 
        % Return Z, Js
    end
end