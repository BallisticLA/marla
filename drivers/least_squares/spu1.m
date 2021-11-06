function[y, log] = spu1(A, c, sampling_factor, tol, iter_lim, logging)
%{
    SVD-based sketch-and-precondition for underdetermined least squares
        min ||y||
        s.t. A' y = c
    where A is tall.
%}
        if logging == 0
            %disp('Optional parameter for logging detailed information has not been passed.'); 
        end

        [n_rows, n_cols] = size(A);
        d = lstsq_dim_checks(sampling_factor, n_rows, n_cols);

        % Sketch the data matrix.
        if logging, tic, end
        % By default, a SJLT sketching matrix is used.
        % Alternative choices are present in '../../Utils/Sketching_Operators'.
        addpath('../../utils/sketching_operators');
        Omega = sjlt(double(d), double(n_rows), double(8));
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
        [x, y, resid_vec] = pcg(A, zeros(n_rows, 1), c, 0.0,...
            tol, iter_lim, M, zeros(n_cols,1));
        if logging, log.t_iterate = toc; end 

        if logging, log.errors = resid_vec; end
end