function SA = srht(A, d, s)
    %Outputs the sketch S*A where S is the subsampled randomized discrete
    %Walsh-Hadamar transform.
    %
    %Utilizes Matlab's built-in fwht() function.
    %
    %Serves an illustrative matter, rather than an efficient
    %implementation.
    %
    %Parameters
    %----------
    %A : matrix
    %    matrix to be sketched
    %d : int
    %    target embedding dimension
    %s : int or RandomStream
    %    Controls random number generation
    %Returns
    %-------
    %SA : Sketch of A of size (d, size(A,2))

    s = MarlaRandStream(s);
    m = size(A, 1);

    % Randomly permute rows of A
    p = randperm(s, m);
    SA = A(p, :);

    % Randomly change signs of the rows of A
    sgn = (rand(s, m, 1) < .5) * 2 - 1;
    SA = sgn .* SA;
    
    % WHT are only defined for even dimensions
    m = 2^(ceil(log2(m)));

    % Applying an m-step WHT
    SA = sqrt(m)*fwht(SA, m);
    
    % Random subsampling of the transformed output
    idx = sort(randsample(s, m, d));
    SA = SA(idx, :);

    % Multiplying by a constanat
    SA = SA * sqrt(m / d);
end
