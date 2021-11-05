function[res, log] = spu1(A, b, sampling_factor, tol, iter_lim, logging)
%{
    SVD-based sketch-and-precondition for underdetermined least squares
        min ||y||
        s.t. A' y = c
    where A is tall.
    #TODO: write proper docstring
%}
        if logging == 0
            %disp('Optional parameter for logging detailed information has not been passed.'); 
        end

        [n_rows, n_cols] = size(A);
        d = dim_checks(sampling_factor, n_rows, n_cols);

        % Sketch the data matrix.
        if logging, tic, end
        % By default, a SJLT sketching matrix is used.
        % Alternative choices are present in '../../Utils/Sketching_Operators'.
        addpath('../../utils/sketching_operators');
        Omega = sjlt(d, num_rows, 8);
        A_ske = Omega * A;
        if logging, log.t_sketch = toc; end

        % Factor the sketch. 
        %  We also measure the time to scale the right singular vectors
        %   as needed for the preconditioner.
        addpath('../../comps/preconditioning/');
        if logging, tic, end
        [M, ~, ~, ~] = svd_right_precond(A_ske);
        if logging, log.t_factor = toc; end

        % Iterative phase
        if logging, tic, end
        addpath('../../comps/itersaddle/');
        % Iterative solver - needs implementation.
        %res = pcss2(A, b, None, c, 0.0, tol, iter_lim, M, False, None);
        if logging, log.t_iterate = toc; end 

        % Need additional 'log' data.
end

% Helper routine. 
function [d] = dim_checks(sampling_factor, num_rows, num_cols)
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