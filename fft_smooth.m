function smooth = fft_smooth(data, threshold)
    n = length(data);
    cutoff = floor(n*(1 - threshold)/2);
    f = fft(data);
    f(n/2+1-cutoff:n/2+cutoff) = zeros(2*cutoff,1);
    smooth = real(ifft(f));
end