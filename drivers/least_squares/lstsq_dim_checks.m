function [d] = lstsq_dim_checks(sampling_factor, num_rows, num_cols)
%{
    Small helper routine, used in lsqr algorithms.
    Takes in integer values of the size of input data and
    a floating point sampling factor, rturns the integer leading dimension
    of the sketching operator. 
%}
    assert(num_rows >= num_cols);
    d = cast((sampling_factor * num_cols), 'uint32');
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
