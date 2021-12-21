function[Out1, Out2, log] = osid1(A, k, over, p, axis, s, logging)
%{
    Sketch + QRCP approach to ID
    See Voronin & Martinsson, 2016, Section 5.1.

    Column ID returns an approximation of the form:
    A ~= (A(:, Out1(1 : k))) * Out2

    Row ID returns an approximation of the form:
    A ~= Out2 * (A(Out1(1 : k), :))

    Important note: 
    When calling a routine, use:
    addpath('../../comps/interpolative/');
    addpath('../../comps/rangefinders/');
%}
    s = MarlaRandStream(s);
    % Relies on a computational routine qrcp_osid
    if axis == 0
        % Row ID
        [S, log] = rs1(A, k + over, p, s, logging);
        Y = A * S;
        [Out1, Out2] = qrcp_osid(Y, k, 0);
    elseif axis == 1
        % Column ID
        [S, log] = rs1(A', k + over, p, s, logging);
        Y = S' * A;
        [Out1, Out2] = qrcp_osid(Y, k, 1);
    end
end