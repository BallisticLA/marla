function [x_ske] = sso1(A, b, sampling_factor, seed) 
%{
    A sketch-and-solve approach to overdetermined ordinary least squares.
    Uses direct method to solve the overdetermined problem.
    ----------
    The sketch-and-solve approach is attributed to a 2006 paper by Sarlos:
    "Improved approximation algorithms for large matrices via random
    projections." An introduction and summary of this approach can be found
    in [MT:2020, Sections 10.2 -- 10.3].

    Important note: 
    ---------------
    Before calling this routine, use:
    addpath('../../utils') - for MatrlaRandStream.m
    addpath('../../Utils/Sketching_Operators')
%}
    seed = MarlaRandStream(seed);
    [num_rows, num_cols] = size(A);
    d = lstsq_dim_checks(sampling_factor, num_rows, num_cols);
    % By default, a SJLT sketching matrix is used.
    % Alternative choices are present in '../../Utils/Sketching_Operators'.
    Omega = sjlt(d, num_rows, 8, seed);
    A_ske = Omega * A;
    b_ske = Omega * b;
    % Solving A_ske * x_ske = b_ske
    [Q, R] = qr(A_ske, 0);
    x_ske = R \ (Q' * b_ske);
end
