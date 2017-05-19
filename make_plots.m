% % make nice plots


%% single and double pulse plots

% load 100F;
load 100C;

% plot of first four double pulses

figure;
plot(xout, ud_out(1:end-1, :));
ylabel('u(x)'); xlabel('x');
legend('2(2)', '2(3)', '2(4)', '2(5)');

% plot of single pulse
figure;
plot(xout, uout(1:end-1));
ylabel('u(x)'); xlabel('x');

%% phase portrait plots

load 100C_timestep;

figure;
hold on;

% unstable point (saddle) and stable point (center)
scatter(6.8295,  0, 'X');
scatter(9.5642,  0, '.');
scatter(12.2880, 0, 'X');
scatter(14.9031, 0, '.');

plot(dist2_2(5:end), deriv2_2(5:end));
plot(dist2_4(5:end), deriv2_4(5:end));
plot(dist2_6(5:end), deriv2_6(5:end));
plot(dist2_8(5:end), deriv2_8(5:end));

plot(dist1_6(5:end), deriv1_6(5:end));

plot(dist3_2(5:end), deriv3_2(5:end));
plot(dist3_4(5:end), deriv3_4(5:end));

plot(dist4_2(5:end), deriv4_2(5:end));
plot(dist4_4(5:end), deriv4_4(5:end));
plot(dist4_6(5:end), deriv4_6(5:end));
plot(dist(5:end), deriv(5:end))


%% exp decay plots

load 100F;
index = 1;

x = xout;
u = uout;
tail_fn    = eigenfunctions(:,(index*2) - 1);
plot_start = length(x)/2 + 35;
plot_end   = length(x) - 20;
xplot = x(plot_start:plot_end);
yplot = log( abs(tail_fn(plot_start:plot_end)) );
c = u(end);

% % roots of linearization about zero solution
% % use for pulses
% nu = roots([1 0 -1 0 par.c]);
% decay = abs(real(nu(1)));

% roots of eigenvalue problem at lambda
% use for eigenfunctions
lambda = l5 * 1i;
nu = roots([-1 0 1 0 -c lambda ]);
decay = abs(max( real( nu(find(real(nu) < -1e-10)))));

plot_name = ['Log of eigenfunction, Double Pulse 2'];
plot_params = ['  c = ', num2str(c),'  lambda = ', num2str(lambda)];
plot_config = [config.method, '  N = ',num2str(length(xout))];
plot_title = [plot_name, plot_params];

marker_size = 5;
scatter(xplot, tail_fn(plot_start:plot_end).*exp(decay*xplot), marker_size, 'filled');



