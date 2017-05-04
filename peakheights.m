% plots heights of peaks of single pulses as a function of c

load ucKdV_Fourier_256;

speeds = [];
heights = [];
for index = 1:length(uc)
    speeds = [speeds; uc(end, index)];
    [pks, locs] = findpeaks( uc(1:end-1, index));
    heights = [heights; max(pks)];
end

bestfit = fit(speeds, heights, 'poly1');
slope =  num2str(abs(bestfit.p1));

plot(speeds, heights);
xlabel('speed (c)');
ylabel('peak height');
title(strcat('Peak height versus speed for KdV5: slope ',slope));