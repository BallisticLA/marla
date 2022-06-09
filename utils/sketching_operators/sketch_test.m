% Script to test correctness of the sketches by verifying that they provide
% subspace embeddings of increasing quality as the embedding dimension is
% increased.

% Sketch to test
sketch = "sjlt"; % Options: "sjlt", "srdft", "srht", "srdct" 

% Generate test matrices
m = 10000;
n = 10;
no_test_vecs = 100;
A = randn(m, n);
X = randn(n, no_test_vecs);

% Compute column norms of A*X before sketching
col_norms = vecnorm(A*X);

% Embedding dimensions to test
embedding_dims = [20, 100, 500, 1000, 2000, 4000, 8000, 10000];

% Compute maximum distortion when sketching
random_seed = 1;
max_distortion = zeros(size(embedding_dims));
for k = 1:length(embedding_dims)
    % Use SJLT
    if sketch == "sjlt"
        if k == 1
            disp("Using sparse JLT sketch")
        end
        Omega = sjlt(embedding_dims(k), m, 8, random_seed);
        SA = Omega * A;
    
    % Use SRDFT
    elseif sketch == "srdft"
        if k == 1
            disp("Using subsampled randomized FFT sketch")
        end
        SA = srdft(A, embedding_dims(k), random_seed);
    
    % Use SRDCT
    elseif sketch == "srdct"
        if k == 1
            disp("Using subsampled randomized DCT sketch")
        end
        SA = srdct(A, embedding_dims(k), random_seed);
    
    % Use SRHT
    elseif sketch == "srht"
        if k == 1
            disp("Using subsampled randomized Hadamard transform sketch")
        end
        SA = srht(A, embedding_dims(k), random_seed);

    % Handle incorrect sketch
    else
        error("Invalid sketch.")
    end

    col_norms_sketched = vecnorm(SA*X);
    max_distortion(k) = max(abs(col_norms_sketched./col_norms - 1));
end

% Plot results
figure
plot(embedding_dims, max_distortion, 'LineWidth', 2, 'Marker', 'o')
xlabel('Embedding dimension')
ylabel('Max distortion')
title("Test of subspace embedding property for " + upper(sketch))
box_dim = [.4 .7 .0 .0];
str = "Distortion should approach zero as \newlineembedding " + ...
    "dimension increases";
annotation('textbox', box_dim,'String', str,'FitBoxToText','on');
shg
