%{
A sketch-and-solve approach to overdetermined ordinary least squares.
----------
The sketch-and-solve approach is attributed to a 2006 paper by Sarlos:
"Improved approximation algorithms for large matrices via random
projections." An introduction and summary of this approach can be found
in [MT:2020, Sections 10.2 -- 10.3].
%}
function [x_ske] = lsqr_sas1(A, b, sampling_factor) 
    [num_rows, num_cols] = size(A);
    d = dim_checks(sampling_factor, num_rows, num_cols);
    Omega = randn(d, num_rows);
    A_ske = Omega * A;
    b_ske = Omega * b;
    x_ske = lsqr(A_ske, b_ske);
end


function [d] = dim_checks(sampling_factor, num_rows, num_cols)
    assert(num_rows >= num_cols);
    d = cast((sampling_factor * num_cols), 'uint16');
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