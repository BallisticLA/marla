function [x_ske] = sso1(A, b, sampling_factor) 
%{
A sketch-and-solve approach to overdetermined ordinary least squares.
Uses direct method to solve the overdetermined problem.
----------
The sketch-and-solve approach is attributed to a 2006 paper by Sarlos:
"Improved approximation algorithms for large matrices via random
projections." An introduction and summary of this approach can be found
in [MT:2020, Sections 10.2 -- 10.3].
%}
    [num_rows, num_cols] = size(A);
    d = dim_checks(sampling_factor, num_rows, num_cols);
    % By default, a SJLT sketching matrix is used.
    % Alternative choices are present in '../../Utils/Sketching_Operators'.
    addpath('../../Utils/Sketching_Operators')
    Omega = sjlt(d, num_rows, 8);
    A_ske = Omega * A;
    b_ske = Omega * b;
    % Solving A_ske * x_ske = b_ske
    [Q, R] = qr(A_ske, 0);
    x_ske = R \ (Q' * b_ske);
end

% Helper routine. 
function [d] = dim_checks(sampling_factor, num_rows, num_cols)
    assert(num_rows >= num_cols);
    d = cast((sampling_factor * num_cols), 'uint8');
    if d > num_rows
        fprintf(['The embedding dimension "d" should not be larger than the', ... 
        'number of rows of the data matrix. Here, an embedding dimension', ...
        'of d=%d has been requested for a matrix with only %d rows.', ...
        'We will proceed by setting d=%d. This parameter choice will', ...
        'result in a very inefficient algorithm!'], d, num_rows, num_rows);
        d = num_rows; 
    end
    assert(d >= num_cols);
end