% One-sided rank-k ID of Y using QR with column pivoting
function [Out1, Out2] = qrcp_osid(Y, k, axis) 
    if axis == 0 
        % Row ID
        [m, ~] = size(Y);
        [~, R, Js] = qr(Y', 'vector');
        I = eye(m, m);
        Z = I(:, Js) * [eye(k, k) R(:, 1 : k) \ R(:, k + 1 : m)]';
        % Returning Z, Js 
        Out1 = Z;
        Out2 = Js;
    elseif axis == 1
        % Column ID
        [Is, Z] = qrcp_osid(Y', k, 0);
        % Returning X, Is
        Out1 = Z';
        Out2 = Is;
    end
end