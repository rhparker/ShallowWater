% semilog plot

N = [10.9375 17.3828 23.6528 29.8828]'
vals = [0.0293 0.0019 0.0006 0.0003]'
xplot = log(N);
yplot = log(vals);

figure;
hold on;
bestfit = fit(xplot, yplot, 'poly1');
scatter(xplot, yplot, 45, 'filled');
plot(xplot, bestfit(xplot));
title(strcat('Log imaginary part of eigenvalue vs spacing of peaks in double pulse, slope = ',num2str(bestfit.p1)));
xlabel('peak spacing');
ylabel('im(eigenvalue)');