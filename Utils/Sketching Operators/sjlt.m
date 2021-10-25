function[Omega] = sjlt(num_rows, num_cols, nnz)
    % Default choice for number of nonzero entries in each column. 
    % nnz = 8;
    
    if num_rows > num_cols
        rows = double.empty;
        nnz = min(num_cols, nnz);
        bad_size = num_rows < nnz;
        if bad_size
            fprintf(['Cannot set %d nonzeros per column for columns ', ...
            'of length %d. Sampling indices with replacement instead.\n\n'], nnz, num_rows);
        end
        
        for i = 1: num_cols
            row = datasample(1:num_rows, nnz, 'Replace', bad_size);
            rows = cat(2, rows, row);
        end
        cols = repelem(1:num_cols, nnz);
        
        % Values for each row and col.
        vals = ones(1, num_cols * nnz);
        vals(rand(1, num_cols * nnz) <= 0.5) = -1;
        vals = vals / sqrt(nnz);   
        % Wrap up.
        Omega = sparse(rows, cols, vals);
    else
        Omega = sjlt(num_cols, num_rows);
        Omega = Omega';
    end
end