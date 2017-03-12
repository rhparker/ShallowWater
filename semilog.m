% semilog plot

pkdist = [6.25 10.1562 13.6719]'
% pkdist = [1 2 3]'
vals   = [0.6423 0.0215 0.0010]'
xplot = log(pkdist);
yplot = vals;

figure;
hold on;
bestfit = fit(xplot, yplot, 'poly1');
scatter(xplot, yplot, 45, 'filled');
plot(xplot, bestfit(xplot));
title(strcat('Log imaginary part of eigenvalue vs spacing of peaks in double pulse, slope = ',num2str(bestfit.p1)));
xlabel('peak spacing');
ylabel('im(eigenvalue)');