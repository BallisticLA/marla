function[Omega] = sjlt(num_rows, num_cols, nnz, s)
    %Generates a sketching operator utilizing Sparse Johnson-Lindenstrauss
    %transform.
    %
    %Parameters
    %----------
    %num_rows : int
    %    number of rows of embedding operator
    %num_cols : int
    %    number of columns of embedding operator
    %nnz : int
    %    number of nonzeros in each column (if num_cols > num_rows) or each 
    %    row (if num_rows >= num_cols)
    %s : int or RandomStream
    %    Controls random number generation
    %Returns
    %-------
    %Omega : Matlab sparse matrix

    % Default choice for number of nonzero entries in each column. 
    % nnz = 8;

    s = MarlaRandStream(s);
    % Ensuring type compatabilities. 
    num_cols = cast(num_cols, 'double');
    num_rows = cast(num_rows, 'double');
    nnz = cast(nnz, 'double');
    
    if num_cols >= num_rows
        rows = double.empty;
        nnz = min(num_cols, nnz);
        bad_size = num_rows < nnz;
        if bad_size
            fprintf(['Cannot set %d nonzeros per column for columns ', ...
            'of length %d. Sampling indices with replacement instead.\n\n'], nnz, num_rows);
        end
        
        for i = 1: num_cols
            row = datasample(s, 1:num_rows, nnz, 'Replace', bad_size);
            rows = cat(2, rows, row);
        end
        cols = repelem(1:num_cols, nnz);
        % Values for each row and col.
        vals = ones(1, num_cols * nnz);
        vals(rand(s, 1, num_cols * nnz) <= 0.5) = -1;
        vals = vals / sqrt(cast(nnz, 'double'));
        % Wrap up.
        Omega = sparse(rows, cols, vals);
    else
        Omega = sjlt(num_cols, num_rows, nnz, s);
        Omega = Omega';
    end
end
