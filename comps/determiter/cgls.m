function[x, iter, resid_vec] = cgls(A, b, tol, iter_lim, precond, x)
    r = b - A * x;
    b_nrm = norm(b);
    p = precond' * (A' * r);
    s = p;
    curr_tol = norm(s, 2)^2;
    stop_tol = (tol*b_nrm)^2;
    resid_vec = [curr_tol];
    iter = 0;
    while curr_tol > stop_tol
        t = precond * p;
        q = A * t;
        alpha = curr_tol / norm(q, 2)^2;
        x = x + alpha * t;
        r = r - alpha * q;
        s = precond' * (A' * r);
        prev_tol = curr_tol;
        resid_vec = [resid_vec, prev_tol]; %#ok<AGROW>
        curr_tol = norm(s, 2)^2;
        beta = curr_tol / prev_tol;
        p = s + beta * p;
        iter = iter + 1;
        if iter == iter_lim
            break
        end
    end
end