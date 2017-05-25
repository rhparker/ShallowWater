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

%% eigenfunction plots

load 100F;
% load 100C;

% % double pulse 2(2) or 2(4)
figure;
% index = 1;
index = 3;
plot(xout, eigenfunctions(:,2*index-1 : 2*index));
ylabel('v(x)'); xlabel('x');
legend(['lambda = ',num2str(eigenvalues(1,index))],['lambda = ',num2str(eigenvalues(2,index))])

% double pulse 2(3) or 2(5)
% index = 2;
index = 4;
% before rotation
% real part
figure;
plot(xout, real(eigenfunctions(:,2*index-1 : 2*index)));
ylabel('Re v(x)'); xlabel('x');
legend(['lambda = ',num2str(eigenvalues(1,index))],['lambda = ',num2str(eigenvalues(2,index))])

% imaginary part
figure;
plot(xout, imag(eigenfunctions(:,2*index-1 : 2*index)));
ylabel('Im v(x)'); xlabel('x');
legend(['lambda = ',num2str(eigenvalues(1,index))],['lambda = ',num2str(eigenvalues(2,index))])

% after rotation
% real part
figure;
plot(xout, real(eigenfunctions_after(:,index/2)));
ylabel('Re v(x)'); xlabel('x');
legend(['lambda = ',num2str(eigenvalues_after(index/2))]);

% imaginary part
figure;
plot(xout, imag(eigenfunctions_after(:,index/2)));
ylabel('Im v(x)'); xlabel('x');
legend(['lambda = ',num2str(eigenvalues_after(index/2))]);

%% eigenvalue plots, Fourier

% Fourier
load 100F;
c = uout(end);

% pulse 2(2)
figure;
plot(lambda1, '.');
axis([-0.6 0.6 -1e6 1e6]);

% pulse 2(3)
figure;
plot(lambda2, '.');
axis([-1e-7 1e-7 -1e6 1e6]);

% pulse 2(3) essential spectrum in exp weighted space
figure;
hold on
plot(lambda2_01, '.');
k = linspace(-16, 16, 10001);
ikalpha = 1i * k - alpha;
ess_spec = ikalpha.^5 - ikalpha.^3 + c*ikalpha;
plot(ess_spec);
legend('Spectrum computed with eig', 'Theoretical essential spectrum');

% pulse 2(3) zoom of exp weighted eigenvalues
figure;
plot(lambda2_01, '.');
axis([-2 0.1 -3 3]);

%% eigenvalue plots, Chebyshev

% Fourier
load 100C;
c = uout(end);

% pulse 2(2), shows absolute spectrum
figure;
plot(lambda1, '.');
axis([-1e4 0.6 -1e5 1e5]);

% pulse 2(2), zoom on eigenvalues
figure;
plot(lambda1, '.');
axis([-0.6 0.6 -1e6 1e6]);

% pulse 2(3)
% figure;
% plot(lambda2, '.');
% axis([-1e-7 1e-7 -1e6 1e6]);


%% phase portrait plots

load 100C_timestep;

figure;
hold on;

% unstable point (saddle) and stable point (center)
scatter(6.8295,  0, 'X');
scatter(9.5642,  0, '.');
scatter(12.2880, 0, 'X');
scatter(14.9031, 0, '.');

% plot(dist1_6(3:400), deriv1_6(3:400));
% plot(dist1_090(3:end), deriv1_090(3:end));
% plot(dist1_095(3:end), deriv1_095(3:end));
% plot(dist1_101(3:end), deriv1_101(3:end));
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

% plot(dist3_2(5:end), deriv3_2(5:end));
% plot(dist3_4(5:end), deriv3_4(5:end));
% 
% plot(dist4_2(5:end), deriv4_2(5:end));
% plot(dist4_4(5:end), deriv4_4(5:end));
% plot(dist4_6(5:end), deriv4_6(5:end));

xlabel('peak distance');
ylabel('derivative of peak distance');


%% eigenvalue decay plots

load 100F;
c = uout(end);

% roots of linearization about 0 solution
nu = roots([1 0 -1 0 c]);
decay = abs(real(nu(1)));
spacing = abs(imag(nu(1)));

% eigenvalues of integrated (4th order) operator
% look at eigenvalue near 0
figure;
marker_size = 15;
hold on;
xplot = spacing*[0 1 2 3]';
yplot = log(abs(int_evals));
scatter(xplot, yplot, marker_size);
% best fit line
bestfit_int = fit(xplot, yplot, 'poly1');
plot(xplot, bestfit_int(xplot));
xlabel('x'); 
ylabel('log(v)');

% eigenvalues of 5th order operator
% here need to look at abs value
figure;
marker_size = 15;
hold on;
xplot = spacing*[0 1 2 3]';
yplot = log(abs(eigenvalues(1,:))');
scatter(xplot, yplot, marker_size);
% best fit line
bestfit = fit(xplot, yplot, 'poly1');
plot(xplot, bestfit(xplot));
xlabel('x'); 
ylabel('log(v)');


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




