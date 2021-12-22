function[Z, Is, X, Js, log] = tsid1(A, k, over, p, s, logging)
%{
    Computes double-sided Interpolative Decomposition of matrix A.
    Rank-k approximation of A is then present as:
    A ~= Z * (A(Is(1 : k), Js(1 : k))) * X.

    Obtain a one-sided ID by any means, then deterministically extend
    to a two-sided ID.
    Using OSID1 would make this a "Sketch + QRCP" approach to double ID,
    as described in Voronin & Martinsson, 2016, Sections 2.4 and 4.
    
    Important note:
    Before running an algorithm, use:     
    addpath('../../comps/interpolative/');
%}
    s = MarlaRandStream(s);
    % Relies on a computational routine qrcp_osid
    if size(A, 1) > size(A, 2)
        [X, Js, log] = osid1(A, k, over, p, 1, s, logging);
        [Z, Is] = qrcp_osid(A(:, Js), k, 0);
    else
        [Z, Is, log] = osid1(A, k, over, p, 0, s, logging);
        [X, Js] = qrcp_osid(A(Is, :), k, 1);
    end
end