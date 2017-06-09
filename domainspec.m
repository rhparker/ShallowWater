% spectrum of a whole bunch of things with domain length varied

% load 100F_domain_768;
% 
% % don't want symmetry enforcement for eigenvalues
% config_nosymm = config;
% config_nosymm.symmetry = 'none';
% config_nosymm.form = 'nonintegrated';
% 
% output_lambda = [];
% total = length(output_L);
% 
% % find spectrum for all the domain lengths
% for index = 1:total
%     x = output_x(:, index);
%     u = output_u(:, index);
%     L = output_L(index);
%     par.c = u(end);              % wave speed
%     N = length(x);               % current grid size
%     
%     % find eigenvalues
%     % no exponential weight
%     a = 0;
%     [lambda, V, J] = eig_linear(x, u(1:end-1), par, config_nosymm, 'nonintegrated', a);
%     output_lambda = [ output_lambda lambda];
% end;

target = 0.0691;
cutoff = 0.001;
output_eigens = [];

figure('DefaultAxesFontSize',16);
hold on;
for index = 1:total
    count = length(output_lambda(:,1));
    plot( output_L(index)*ones(1, count), imag(output_lambda(:, index)), '.', 'color', 'blue');
    % find above eigenvalues in unweighted spectrum
    lambda = output_lambda(:,index);
    eigens = lambda(find( abs( abs(imag(lambda)) - target) < cutoff));
    output_eigens = [output_eigens eigens];
    plot( output_L(index)*[1 1], imag(eigens), '.', 'Color', 'red', 'MarkerSize', 10);
    % predicted lowest point of essential spectrum
    c = 10;
    essbottom = [pi*c / output_L(index) -pi*c / output_L(index)];
    plot( output_L(index)*[1 1], essbottom, '.', 'Color', 'green', 'MarkerSize', 10);
end
ybound = 1;
axis([output_L(1) output_L(end) -ybound ybound]);
xlabel('domain length (L)');
ylabel('imaginary part of spectrum');

figure('DefaultAxesFontSize',16);
hold on;
xplot = output_L;
yplot = log( abs( real(output_eigens(1,:))));
plot(xplot, yplot, '.', 'MarkerSize', 10);
bestfit = fit(xplot', yplot', 'poly1');
plot(xplot, bestfit(xplot));
xlabel('Domain length (L)');
ylabel('Log of real part of nonzero eigenvalue');
legend('log(real(lambda))',['slope = ',num2str(bestfit.p1)]);

