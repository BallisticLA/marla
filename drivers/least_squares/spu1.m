function[y, log] = spu1(A, c, sampling_factor, tol, iter_lim, logging, seed)
%{
    SVD-based sketch-and-precondition for underdetermined least squares
        min ||y||
        s.t. A' y = c
    where A is tall.

    Important note: 
    ---------------
    Before calling this routine, use:
    addpath('../../utils/sketching_operators');
    addpath('../../comps/preconditioning/');
    addpath('../../comps/itersaddle/');
%}
        if logging.depth == 0 || logging.span == 0
            log_present = 0;
            log.status = 'Optional parameter for logging detailed information has not been passed.'; 
        else
            log_present = 1;
            logging.depth = logging.depth - 1;
        end

        seed = MarlaRandStream(seed);

        [n_rows, n_cols] = size(A);
        d = lstsq_dim_checks(sampling_factor, n_rows, n_cols);

        % Sketch the data matrix.
        if log_present, tic, end
        % By default, a SJLT sketching matrix is used.
        % Alternative choices are present in '../../Utils/Sketching_Operators'.
        Omega = sjlt(double(d), double(n_rows), double(8), seed);
        A_ske = Omega * A;
        if log_present, log.t_sketch = toc; end

        % Factor the sketch. 
        %  We also measure the time to scale the right singular vectors
        %   as needed for the preconditioner.
        if log_present, tic, end
        [M, ~, ~, ~] = svd_right_precond(A_ske);
        if log_present, log.t_factor = toc; end

        % Iterative phase
        if log_present, tic, end
        % Iterative solver - needs implementation.
        [~, y, resid_vec] = pcg(A, zeros(n_rows, 1), c, 0.0,...
            tol, iter_lim, M, zeros(n_cols,1));
        if log_present, log.t_iterate = toc; end 

        if log_present, log.errors = resid_vec; end
end