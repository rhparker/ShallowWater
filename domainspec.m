% spectrum of a whole bunch of things with domain length varied

load 100F_domain;

% don't want symmetry enforcement for eigenvalues
config_nosymm = config;
config_nosymm.symmetry = 'none';
config_nosymm.form = 'nonintegrated';

output_lambda = [];
total = length(output_L);

% find spectrum for all the domain lengths
for index = 1:total
    x = output_x(:, index);
    u = output_u(:, index);
    L = output_L(index);
    par.c = u(end);              % wave speed
    N = length(x);               % current grid size
    
    % find eigenvalues
    % no exponential weight
    a = 0;
    [lambda, V, J] = eig_linear(x, u(1:end-1), par, config_nosymm, 'nonintegrated', a);
    output_lambda = [ output_lambda lambda];
end;

figure;
hold on;
for index = 1:total
    plot( index*ones(1, length(x)), imag(output_lambda(:, index)), '.');
end
ybound = 2;
axis([0 total+1 -ybound ybound]);
