% N  = [1000 2000 3000 4000 5000]';
% nu = [.0088 .0022 .0010 .0005 .0003]';
% L  = 50;
% h  = 2*L./N;

% cutoff = length(leigs02);
% index  = 2;
% true_value = 0;
% grids  = leigs02(1:cutoff,1);
% vals   = abs( real( leigs02(1:cutoff,index) ) - true_value);
% 
% grids = integsAvg(:,1);
% vals  = abs(real(integsAvg(:,4)));
% 
% % diffs = true;
% diffs = false;
% if diffs
%     vals  = abs( vals(2:end) - vals(1:end-1) );
%     grids = grids(1:end-1);
% end
% 
% xplot = log( 2*L./ grids );
% % xplot = grids;
% yplot = log(vals);

%
% tail plots (look for exp decay)
%
load 5double1;
x          = xout;
tail_fn    = ud_out;
plot_start = length(x)/2 + 35;
plot_end   = length(x)-25;
xplot = x(plot_start:plot_end);
yplot = log( abs(tail_fn(plot_start:plot_end)) );
par.c = uout(end);
nu = roots([1 0 -1 0 par.c]);
decay = abs(real(nu(1)));
plot_name = ['Log of single pulse'];
plot_params = ['  c = ', num2str(par.c),'  mu = ',num2str(decay),'   ',];
plot_config = [config.method, '  N = ',num2str(length(xout))];
plot_title = [plot_name, plot_params, plot_config];


% % maxflip plots
% mf = maxflips;
% yplot = log( mf(:,3) );
% xplot = log( mf(:,2) );


figure;
hold on;
marker_size = 10;

% titles and labels
title(plot_title);
ylabel('log of pulse');
xlabel('x');

scatter(xplot, yplot, marker_size, 'filled');

% best fit line 
bestfit = fit(xplot, yplot, 'poly1');
plot(xplot, bestfit(xplot));
title({ plot_title, strcat('slope =  ',num2str(abs(bestfit.p1))) } );

