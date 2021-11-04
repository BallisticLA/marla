function[Out1, Out2] = osid1(A, k, s, p, axis)
%{
    Sketch + QRCP approach to ID
    See Voronin & Martinsson, 2016, Section 5.1.

    Column ID returns an approximation of the form:
    A ~= (A(:, Out1(1 : k))) * Out2

    Row ID returns an approximation of the form:
    A ~= Out2 * (A(Out1(1 : k), :))
%}
    % Relies on a computational routine qrcp_osid
    % Uses rf1 as a default rangefinder, but alternatives are available. 
    addpath('../../comps/interpolative/');
    addpath('../../comps/rangefinders/');
    if axis == 0
        % Row ID
        Y = rf1(A, k + s, p);
        % return X, Is
        [Out1, Out2] = qrcp_osid(Y, k, 0);
    elseif axis == 1
        % Column ID
        Y = rf1(A', k+s, p);
        % return X, Is
        [Out1, Out2] = qrcp_osid(Y, k, 1);
    end
end