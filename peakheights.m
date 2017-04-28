% plots heights of peaks of single pulses as a function of c

load ucKdV_Fourier_256;

speeds = [];
heights = [];
for index = 1:length(uc)
    speeds = [speeds; uc(end, index)];
    [pks, locs] = findpeaks( uc(1:end-1, index));
    heights = [heights; max(pks)];
end

plot(speeds, heights)