function[Out1, Out2, log] = osid2(A, k, over, p, axis, s, logging)
%{
    Sketch + (QRCP skeleton) + (least squares) approach to ID
    See Dong & Martinsson, 2021.

    Column ID returns an approximation of the form:
    A ~= (A(:, Out1(1 : k))) * Out2

    Row ID returns an approximation of the form:
    A ~= Out2 * (A(Out1(1 : k), :))

    Important note: 
    When calling a routine, use:
    addpath('../../comps/rangefinders/');
%}
    s = MarlaRandStream(s);
    if axis == 0
        % Row ID
        [S, log] = rs1(A, k + over, p, s, logging);
        Y = A * S;
        [~, ~, I] = qr(Y', 'vector');
        Out2 = I(1 : k);
        Out2 = Out2(:);
        Out1 = A / A(Out2, :);
    elseif axis == 1
        % Column ID
        [S, log] = rs1(A', k + over, p, s, logging);
        Y = S' * A;
        [~, ~, J] = qr(Y, 'vector');
        Out2 = J(1 : k);
        Out2 = Out2(:);
        Out1 = A(:, Out2) \ A; 
    end
end