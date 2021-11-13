function[Out1, Out2] = osid1(A, k, over, p, axis, s)
%{
    Sketch + QRCP approach to ID
    See Voronin & Martinsson, 2016, Section 5.1.

    Column ID returns an approximation of the form:
    A ~= (A(:, Out1(1 : k))) * Out2

    Row ID returns an approximation of the form:
    A ~= Out2 * (A(Out1(1 : k), :))
%}
    s = MarlaRandStream(s);
    % Relies on a computational routine qrcp_osid
    addpath('../../comps/interpolative/');
    addpath('../../comps/rangefinders/');
    if axis == 0
        % Row ID
        S = rs1(A, k + over, p, s);
        Y = A * S;
        [Out1, Out2] = qrcp_osid(Y, k, 0);
    elseif axis == 1
        % Column ID
        S = rs1(A', k + over, p, s)';
        Y = S * A;
        [Out1, Out2] = qrcp_osid(Y, k, 1);
    end
end