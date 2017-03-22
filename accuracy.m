% N  = [1000 2000 3000 4000 5000]';
% nu = [.0088 .0022 .0010 .0005 .0003]';
% L  = 50;
% h  = 2*L./N;

cutoff = length(leigs02);
index  = 2;
true_value = 0;
grids  = leigs02(1:cutoff,1);
vals   = abs( real( leigs02(1:cutoff,index) ) - true_value);

grids = integsAvg(:,1);
vals  = abs(real(integsAvg(:,4)));

% diffs = true;
diffs = false;
if diffs
    vals  = abs( vals(2:end) - vals(1:end-1) );
    grids = grids(1:end-1);
end

xplot = log( 2*L./ grids );
% xplot = grids;
yplot = log(vals);

% % tail plots 
% plot_start = length(xnew)/2 + 50;
% plot_end   = length(xnew);
% xplot = xnew(plot_start:plot_end);
% yplot = log10( abs(eVecs(plot_start:plot_end)) );

% flip pmaxlots
mf = maxflips;
yplot = log( mf(:,3) );
xplot = log( mf(:,2) );


figure;
hold on;
marker_size = 40;

% titles and labels
plot_title = 'Log max |v(x)| - |v(-x)| vs log L, lambda=0.0215i, Fourier, N=2048';
title(plot_title);
ylabel('log max |v(x)| - |v(-x)|');
xlabel('log L');

scatter(xplot, yplot, marker_size, 'filled');

% best fit line 
bestfit = fit(xplot, yplot, 'poly1');
plot(xplot, bestfit(xplot));
title(strcat(plot_title,', slope = ',num2str(abs(bestfit.p1))));

