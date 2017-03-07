% N  = [1000 2000 3000 4000 5000]';
% nu = [.0088 .0022 .0010 .0005 .0003]';
% L  = 50;
% h  = 2*L./N;

% grids = [5000, 10000, 20000, 40000]';
% diff   = abs( real(leigs02(2:end,1) ))
% change = diff(1:end-1) - diff(2:end);

grids = Neigs02;
diffs  = abs( imag( leigs02(:,1) ) - 0.18  );

xplot = log(200./grids);
yplot = log( diffs );
% yplot = log( abs( leigs(:,1) - 0.605i ) )
figure;
hold on;
bestfit = fit(xplot, yplot, 'poly1');
scatter(xplot, yplot, 45, 'filled');
plot(xplot, bestfit(xplot));
title(strcat('Log abs value of eigenvalue vs Log mesh size (c = 18.6836), order = ',num2str(abs(bestfit.p1))));
xlabel('Log mesh size');
ylabel('Log real part of eigenvalue');

