function w = cheb_weights(N)
    cj = [2; ones(N-2, 1); 2];
    w  = pi ./ (N * cj);
end