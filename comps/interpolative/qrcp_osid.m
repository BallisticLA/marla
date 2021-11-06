% One-sided rank-k ID of Y using QR with column pivoting
function [Out1, Out2] = qrcp_osid(Y, k, axis) 
    if axis == 1
        % n = size(Y,2);
        [~, R, J] = qr(Y,'vector');
        temp = [eye(k, k) R(1:k, 1:k) \ R(1:k, k + 1 : end)];
        X(:, J) = temp;
        Out1 = X;
        Out2 = J(1:k);
        Out2 = Out2(:);
    elseif axis == 0
        [X, Is] = qrcp_osid(Y', k, 1);
        Out1 = X';
        Out2 = Is(:);  % force a column vector
    end
    %if axis == 0 
    %    % Row ID
    %    [m, ~] = size(Y);
    %    [~, R, Js] = qr(Y', 'vector');
    %    Z(:, Js) = [eye(k, k) R(:, 1:k) \ R(:, k + 1 : m)]';
    %    % Returning Z, Js 
    %    Out1 = Z;
    %    Out2 = Js(1:k);
    %elseif axis == 1
    %    % Column ID
    %    [Is, Z] = qrcp_osid(Y', k, 0);
    %    % Returning X, Is
    %    Out1 = Z';
    %    Out2 = Is(1:k);
    %end
end