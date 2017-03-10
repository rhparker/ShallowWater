% N  = [1000 2000 3000 4000 5000]';
% nu = [.0088 .0022 .0010 .0005 .0003]';
% L  = 50;
% h  = 2*L./N;

L = 25;

cutoff = length(leigs02);
index  = 2;
true_value = 0;
grids  = leigs02(1:cutoff,1);
vals   = abs( real( leigs02(1:cutoff,index) ) - true_value);

% diffs = true;
diffs = false;
if diffs
    vals  = abs( vals(2:end) - vals(1:end-1) );
    grids = grids(1:end-1);
end

xplot = log( 2*L./ grids);
yplot = log( vals );

% yplot = log( abs( leigs(:,1) - 0.605i ) )

figure;
hold on;
plot_title = 'log (abs value of) real part of Left Eigenvalue vs log mesh size';
ylabel('log (abs value of) real part of Left Eigenvalue');
bestfit = fit(xplot, yplot, 'poly1');
scatter(xplot, yplot, 45, 'filled');
plot(xplot, bestfit(xplot));
title(strcat(plot_title,', order = ',num2str(abs(bestfit.p1))));
xlabel('Log mesh size');

