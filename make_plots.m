% % make nice plots

%% single pulse plots for different values of c

load 100F;
u100 = uout;
load 75F;
u75 = uout;
load 50F;
u50 = uout;
load 25F;
u25 = uout;

uwaves = [u25 u50 u75 u100];
figure('DefaultAxesFontSize',16);
plot(xout, uwaves(1:end-1,:));
xlabel('x');
ylabel('u(x)');
legend('c = 2.5', 'c = 5.0', 'c = 7.5', 'c = 10');


%% pulse heights vs speed for single pulses

load ucKdV_Fourier_256;
hts = max(uc(1:end-1,:));
speeds = uc(end, :);
figure('DefaultAxesFontSize',16);
plot(speeds, hts);
xlabel('wave speed (c)');
ylabel('peak height of single pulse');
bestfit = fit(speeds', hts', 'poly1');
legend(['slope: ' num2str(bestfit.p1)]);


%% single and double pulse plots

load 100F;

% plot of first four double pulses

ud_four = ud_out(1:end-1, :);
x_four = xout;
figure('DefaultAxesFontSize',16);
plot(xout, ud_four);
ylabel('u(x)'); xlabel('x');
legend('2(1)', '2(2)', '2(3)', '2(4)', '2(5)');

load 100C;

ud_cheb = ud_out(1:end-1, :);
x_cheb = xout;
figure('DefaultAxesFontSize',16);
plot(xout, ud_cheb);
ylabel('u(x)'); xlabel('x');
legend('2(1)', '2(2)', '2(3)', '2(4)', '2(5)');

% do with L = 40
load 100F_40;
ud_four = ud_out(1:end-1, :);
x_four = xout;

load 100C_40;
ud_cheb = udc_out(1:end-1, :);
x_cheb = xcout;

% plot of differences between methods
figure('DefaultAxesFontSize',16);
% interpolate so we can subtract
ud_four_interp = ud_cheb;
ud_cheb_interp = ud_four;
for index = 1:4 
    ud_four_interp(:,index) = spline(x_four, ud_four(:,index), x_cheb);
    ud_cheb_interp(:,index) = spline(x_cheb, ud_cheb(:,index), x_four);
end
plot(x_cheb, ud_cheb - ud_four_interp);
ylabel('Chebyshev u(x) - Fourier u(x)'); xlabel('x');
legend('2(2)', '2(3)', '2(4)', '2(5)');
% plot(x_four, ud_four - ud_cheb_interp);


%%  plot of multipulses

load 100F;
figure('DefaultAxesFontSize',16);
plot(xout, um_1_3(1:end-1));
ylabel('u(x)'); xlabel('x');

figure('DefaultAxesFontSize',16);
plot(xout, um_1_4(1:end-1));
ylabel('u(x)'); xlabel('x');

load 100C;

figure('DefaultAxesFontSize',16);
plot(xout, umc_2_3(1:end-1));
ylabel('u(x)'); xlabel('x');

figure('DefaultAxesFontSize',16);
plot(xout, umc_2_4(1:end-1));
ylabel('u(x)'); xlabel('x');

%% eigenfunction plots

load 100F;
% load 100C;

% % double pulse 2(2) or 2(4)
figure('DefaultAxesFontSize',16);
% index = 1;
index = 3;
plot(xout, eigenfunctions(:,2*index-1 : 2*index));
ylabel('v(x)'); xlabel('x');
legend(['lambda = ',num2str(eigenvalues(1,index))],['lambda = ',num2str(eigenvalues(2,index))])

% double pulse 2(3) or 2(5)
index = 2;
% index = 4;
% before rotation
% real part
figure('DefaultAxesFontSize',16);
plot(xout, real(eigenfunctions(:,2*index-1 : 2*index)));
ylabel('Re v(x)'); xlabel('x');
legend(['Im lambda = ',num2str(imag(eigenvalues(1,index)))],['Im lambda = ',num2str(imag(eigenvalues(2,index)))])

% imaginary part
figure('DefaultAxesFontSize',16);
plot(xout, imag(eigenfunctions(:,2*index-1 : 2*index)));
ylabel('Im v(x)'); xlabel('x');
legend(['Im lambda = ',num2str(imag(eigenvalues(1,index)))],['Im lambda = ',num2str(imag(eigenvalues(2,index)))])

% after rotation
% real part
figure('DefaultAxesFontSize',16);
plot(xout, real(eigenfunctions_after(:,index/2)));
ylabel('Re v(x)'); xlabel('x');
legend(['Im lambda = ',num2str(imag(eigenvalues_after(index/2)))]);

% imaginary part
figure('DefaultAxesFontSize',16);
plot(xout, imag(eigenfunctions_after(:,index/2)));
ylabel('Im v(x)'); xlabel('x');
legend(['Im lambda = ',num2str(imag(eigenvalues_after(index/2)))]);

%% eigenvalue plots, Fourier

% Fourier
load 100F;
c = uout(end);

% pulse 2(2)
figure('DefaultAxesFontSize',16);
plot(lambda1, '.', 'MarkerSize', 10);
axis([-0.6 0.6 -1e6 1e6]);

% pulse 2(3)
figure('DefaultAxesFontSize',16);
plot(lambda2, '.', 'MarkerSize', 10);
axis([-1e-7 1e-7 -1e6 1e6]);

% pulse 2(3) essential spectrum in exp weighted space
figure('DefaultAxesFontSize',16);
hold on
plot(lambda2_01, '.');
k = linspace(-16, 16, 10001);
ikalpha = 1i * k - alpha;
ess_spec = ikalpha.^5 - ikalpha.^3 + c*ikalpha;
plot(ess_spec);
legend('Spectrum computed with eig', 'Predicted essential spectrum');

% pulse 2(3) zoom of exp weighted eigenvalues
figure('DefaultAxesFontSize',16);
plot(lambda2_01, '.', 'MarkerSize', 10);
axis([-2 0.1 -3 3]);

% pulse 3(3,3)
figure('DefaultAxesFontSize',16);
plot(lambda1_triple, '.', 'MarkerSize', 10);
axis([-0.6 0.6 -1e6 1e6]);

% pulse 3(3,3,3)
figure('DefaultAxesFontSize',16);
plot(lambda1_quad, '.', 'MarkerSize', 10);
axis([-0.6 0.6 -1e6 1e6]);

%% eigenvalue plots, Fourier
% "fishing" out eigenvalues from essential spectrum
% after finding in weighed space

load 100F;
% load 250F;

% same y bound for both plots
ybound = 2;

% cutoff for Matlab's find
cutoff = 0.0001;

% pulse 2(3) zoom of exp weighted eigenvalues
figure('DefaultAxesFontSize',16);
hold on;
plot(lambda2_01, '.', 'MarkerSize', 10);
% acutal eigenvalues are on imaginary axis
eigens = lambda2_01( find(abs(real(lambda2_01)) < cutoff ) );
% only want the nonzero eigenvalues
eigens = eigens(find(abs(imag(eigens)) > cutoff ));
plot(eigens, '.', 'MarkerSize', 10, 'Color', 'red');
axis([-2 0.1 -ybound ybound]);

% pulse 2(3)
figure('DefaultAxesFontSize',16);
hold on;
plot(lambda2, '.', 'MarkerSize', 10);
target = abs(imag(eigens(1)));
% find above eigenvalues in unweighted spectrum
eigens_unwt = lambda2(find( abs( abs(imag(lambda2)) - target) < cutoff));
plot(eigens_unwt, '.', 'MarkerSize', 10, 'Color', 'red');
axis([-0.01 0.01 -ybound ybound]);


%% eigenvalue plots, Chebyshev

% Fourier
load 100C;
c = uout(end);

% % pulse 2(2), shows absolute spectrum
% figure('DefaultAxesFontSize',16);
% plot(lambda1, '.');
% axis([-1e4 0.6 -1e5 1e5]);
% 
% % pulse 2(2), zoom on eigenvalues
% figure('DefaultAxesFontSize',16);
% plot(lambda1, '.', 'MarkerSize', 10);
% axis([-0.6 0.6 -1e6 1e6]);

% pulse 2(3)
figure('DefaultAxesFontSize',16);
plot(lambda2, '.');
% axis([-1e-7 1e-7 -1e6 1e6]);

% pulse 2(3), zoom on eigenvalues
figure('DefaultAxesFontSize',16);
plot(lambda2, '.', 'MarkerSize', 10);
axis([-30 0.6 -5 5]);

% % pulse 3(2,2)
% figure('DefaultAxesFontSize',16);
% plot(lambda2_triple, '.', 'MarkerSize', 10);
% axis([-0.1 0.1 -0.2 0.2]);
% 
% % pulse 3(2,2,2)
% figure('DefaultAxesFontSize',16);
% plot(lambda2_quad, '.', 'MarkerSize', 10);
% axis([-0.6 0.6 -0.2 0.2]);

%% eigenvalues vs speed c

speeds = [5 7.5 10 15 20 25]';
logspeeds = log(speeds);
lambda1 = [0.1363 0.2730 0.4352 0.8157 1.2535 1.7360]';
lambda2 = [0.0190 0.0413 0.0691 0.1367 0.2168 0.3067]';

figure('DefaultAxesFontSize',16);
hold on;
plot(log(speeds), [log(lambda1) log(lambda2)], '.', 'MarkerSize',10);
xbestfit2 = fit(logspeeds, log(lambda2), 'poly1');
plot(logspeeds, bestfit1(logspeeds));
plot(logspeeds, bestfit2(logspeeds));
xlabel('log of speed (log c)');
ylabel('log of abs value of eigenvalue (log |lambda|)');
slope1 = ['y = ' num2str(bestfit1.p1) ' x + ' num2str(bestfit1.p2) ];
slope2 = ['y = ' num2str(bestfit2.p1) ' x + ' num2str(bestfit2.p2) ];
legend('double pulse 2(2)', 'double pulse 2(3)', slope1, slope2, 'Location','northwest');


%% decay rate alpha vs c

cvals = linspace(1, 100, 100);
alphavals = zeros(1, length(cvals));
for index = 1:length(cvals);
    nu = roots([1 0 -1 0 cvals(index)]);
    alpha = abs(real(nu(1)));
    alphavals(index) = alpha;
end

% figure('DefaultAxesFontSize',16);
% plot(cvals, alphavals, '.');
% xlabel('speed (c)');
% ylabel('alpha');
% 
% 

figure('DefaultAxesFontSize',16);
logc = log(cvals);
logalpha = log(alphavals);
plot(logc, logalpha, '.');
xlabel('log of speed (log c)');
ylabel('log(alpha)');


figure('DefaultAxesFontSize',16);
hold on;
xplot = logc( 30 : end);
yplot = logalpha( 30 : end);
plot(xplot, yplot, '.');
bestfit = fit(xplot', yplot', 'poly1');
plot(xplot, bestfit(xplot));
xlabel('log of speed (log c)');
ylabel('log(alpha)');
slope = ['y = ' num2str(bestfit.p1) ' x + ' num2str(bestfit.p2) ];
legend('log(alpha)', slope, 'Location','northwest');

%% integral of eigenfunctions

load 100F;
N = length(xout);               % current grid size
L = ceil(abs(xout(1)));         % current domain length
h = 2*L/N;                      % current grid spacing
Fintegs = h * sum(eigenfunctions);
% only need to do this for one eigenvalue per double pulse
Fintegs = Fintegs([1 3 5 7]);

load 100C;



%% phase portrait plots

load 100C_timestep;

figure('DefaultAxesFontSize',16);
hold on;

% unstable point (saddle) and stable point (center)
scatter(6.8295,  0, 'X');
scatter(9.5642,  0, '.');
scatter(12.2880, 0, 'X');
scatter(14.9031, 0, '.');

% remove this block for zoom
plot(dist1_6(3:400), deriv1_6(3:400));
plot(dist1_090(3:end), deriv1_090(3:end));
plot(dist1_095(3:end), deriv1_095(3:end));
plot(dist1_101(3:end), deriv1_101(3:end));

plot(dist1_105(3:150), deriv1_105(3:150));
plot(dist1_110(3:200), deriv1_110(3:200));

plot(dist2_2(5:end), deriv2_2(5:end));
plot(dist2_4(5:end), deriv2_4(5:end));
plot(dist2_6(5:end), deriv2_6(5:end));
plot(dist2_8(5:end), deriv2_8(5:end));
plot(dist2_095(5:end), deriv2_095(5:end));
plot(dist2_092(5:end), deriv2_092(5:end));
plot(dist2_091(5:end), deriv2_091(5:end));
plot(dist2_0915(5:end), deriv2_0915(5:end));

% remove this block for zoom
plot(dist3_2(5:end), deriv3_2(5:end));
plot(dist3_4(5:end), deriv3_4(5:end));
plot(dist4_2(5:end), deriv4_2(5:end));
plot(dist4_4(5:end), deriv4_4(5:end));
plot(dist4_6(5:end), deriv4_6(5:end));

xlabel('peak distance');
ylabel('derivative of peak distance');

%% plots of timestepping at stationary points

load 100C_timestep;

% double pulse 2

figure('DefaultAxesFontSize',16);
plot(time2_0, dist2_0);
xlabel('time');
ylabel('peak distance');

figure('DefaultAxesFontSize',16);
plot(dist2_0, deriv2_0);
xlabel('peak distance');
ylabel('derivative of peak distance');

[pks, locs] = findpeaks(dist2_0);
gaps = locs(2:end) - locs(1:end-1);
mean_gap = mean(gaps);
timestep = time2_0(2) - time2_0(1);
period = mean_gap * timestep;


%% eigenvalue decay plots

load 100F;
marker_size = 15;
c = uout(end);

% roots of linearization about 0 solution
nu = roots([1 0 -1 0 c]);
decay = abs(real(nu(1)));
spacing = abs(imag(nu(1)));

% eigenvalues of integrated (4th order) operator
% look at eigenvalue near 0
figure('DefaultAxesFontSize',16);
% marker_size = 15;
hold on;
xplot = (pi / spacing)*[0 1 2 3 4]';
yplot = log(abs(int_evals));
scatter(xplot, yplot, marker_size);
% best fit line
bestfit_int = fit(xplot, yplot, 'poly1');
plot(xplot, bestfit_int(xplot));
xlabel('distance between pulses (normalized)'); 
ylabel('log(lambda)');
legend('log(lambda)',['slope = ',num2str(bestfit_int.p1)]);

% eigenvalues of 5th order operator
% here need to look at abs value
figure('DefaultAxesFontSize',16);
marker_size = 25;
hold on;
xplot = [0 1 2 3]';
yplot = log(abs(eigenvalues(1,:))');
const = abs(eigenvalues(1,1)) / sqrt(spacing);
predictedeigs = const * sqrt(spacing)*exp(-decay * (xplot) * pi / (2 * spacing));
predictedplot = log(predictedeigs);
scatter(xplot, yplot, marker_size);
scatter(xplot, predictedplot, marker_size, 'x');
% best fit line
bestfit = fit(xplot, yplot, 'poly1');
plot(xplot, bestfit(xplot), 'b');
xlabel('n'); 
ylabel('log(lambda)');
legend('abs value of eigenvalue', 'predicted abs value of eigenvalue');


%% exp decay plots

load 100F;

x = xout;
u = uout;
c = u(end);

index = 1;

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
index = 1;
lambda = eigenvalues(1, index);
nu = roots([-1 0 1 0 -c lambda ]);
decay = abs(max( real( nu(find(real(nu) < -1e-10)))));
tail_fn    = real( eigenfunctions(:,(index*2) - 1) );

%
% plot bounds
%

% % these work for single pulses
% plot_start = length(x)/2   + 4;
% plot_end   = length(x)     - 18;

% % these work for eigenfunctions
% for first eigenfunction
plot_start = length(x)/2   + 23;
plot_end   = length(x)     - 65;
 
% % for second eigenfunction
% plot_start = length(x)/2   + 28;
% plot_end   = length(x)     - 55;

% first plot is log of tail of function vs x
figure('DefaultAxesFontSize',16);
hold on;
xplot = x( plot_start:plot_end );
yplot = log( abs(tail_fn(plot_start:plot_end)) );
plot(xplot, yplot, '.');
% best fit line
bestfit = fit(xplot, yplot, 'poly1');
plot(xplot, bestfit(xplot));
xlabel('x'); 
% ylabel('log(u)');
% legend('log of pulse', ['best fit line, slope =',num2str(bestfit.p1)]);
ylabel('log(v)');
legend('log of eigenfunction', ['best fit line, slope =',num2str(bestfit.p1)]);


figure('DefaultAxesFontSize',16);
marker_size = 15;
yplot_scale = tail_fn(plot_start:plot_end).*exp(decay*xplot);
scatter(xplot, yplot_scale, marker_size, 'filled');
xlabel('x');
% ylabel('u exp(alpha x)');
ylabel('v exp(alpha x)');




