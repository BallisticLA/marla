function[Z, Is, X, Js] = tsid1(A, k, s, p)
%{
    Computes double-sided Interpolative Decomposition of matrix A.
    Rank-k approximation of A is then present as:
    A ~= Z * (A(Is(1 : k), Js(1 : k))) * X.

    Obtain a one-sided ID by any means, then deterministically extend
    to a two-sided ID.
    Using OSID1 would make this a "Sketch + QRCP" approach to double ID,
    as described in Voronin & Martinsson, 2016, Sections 2.4 and 4.
%}
    % Relies on a computational routine qrcp_osid
    addpath('../../comps/interpolative/');
    if size(A, 1) > size(A, 2)
        [X, Js] = osid1(A, k, s, p, 0);
        [Z, Is] = qrcp_osid(A(:, Js), k, 0);
    else
        [Z, Is] = osid(A, k, s, p, 1);
        [X, Js] = qrcp_osid(A(Is, :), k);
    end
end