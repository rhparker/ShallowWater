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

% plot(dist1_6(5:end), deriv1_6(5:end));

% plot(dist3_2(5:end), deriv3_2(5:end));
% plot(dist3_4(5:end), deriv3_4(5:end));
% 
% plot(dist4_2(5:end), deriv4_2(5:end));
% plot(dist4_4(5:end), deriv4_4(5:end));
% plot(dist4_6(5:end), deriv4_6(5:end));

plot(dist2_095(5:end), deriv2_095(5:end));
plot(dist2_092(5:end), deriv2_092(5:end));

plot(dist(5:end), deriv(5:end))


%% exp decay plots

load 100F;

x = xout;
u = uout;
c = u(end);

index = x;

% % % single or double pulses
% % roots of linearization about zero solution
% % use for pulses
% nu = roots([1 0 -1 0 c]);
% decay = abs(real(nu(1)));
% spacing = abs(imag(nu(1)));
% tail_fn = uout;

% % eigenfunction
% roots of eigenvalue problem at lambda
% use for eigenfunctions
index = 2;
lambda = eigenvalues(1, index);
nu = roots([-1 0 1 0 -c lambda ]);
decay = abs(max( real( nu(find(real(nu) < -1e-10)))));
tail_fn    = real( eigenfunctions(:,(index*2) - 1) );

%
% plot bounds
%

% these work for single pulses
plot_start = length(x)/2   + 4;
plot_end   = length(x)     - 18;

% these work for eigenfunctions
% for first eigenfunction
plot_start = length(x)/2   + 23;
plot_end   = length(x)     - 65;

% for second eigenfunction
plot_start = length(x)/2   + 28;
plot_end   = length(x)     - 55;

% first plot is log of tail of function vs x
figure;
hold on;
xplot = x( plot_start:plot_end );
yplot = log( abs(tail_fn(plot_start:plot_end)) );
plot(xplot, yplot, '.');
% best fit line
bestfit = fit(xplot, yplot, 'poly1');
plot(xplot, bestfit(xplot));
xlabel('x'); 
ylabel('log(u)');
% ylabel('log(v)');


figure;
marker_size = 15;
yplot_scale = tail_fn(plot_start:plot_end).*exp(decay*xplot);
scatter(xplot, yplot_scale, marker_size, 'filled');
xlabel('x');
ylabel('u exp(alpha x)');
% ylabel('v exp(alpha x)');




