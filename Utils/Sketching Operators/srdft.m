function[Omega] = srdft(A, d)
    [m, n] = size(A);
    % Generating a random sign vector
    sgn = (rand(1, m) < .5) * 2 - 1;
    % Randomly changing signs of columns of A
    A = bsxfun(@times, A, sgn);
    
    % Applying FFT
    Omega = (fft(A));
    
    % Random subsampling of the transform output
    idx = sort(randsample(m, d));
    Omega = Omega(idx, :);
    % Multiplying by a constanat
    Omega = Omega * (sqrt(m / d));
    % Optional random column permutation
    idx = randperm(n);
    Omega = Omega(:,idx);
end