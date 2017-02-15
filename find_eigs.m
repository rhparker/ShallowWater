%% eigenvalues of integrated operator

% % single pulse
% load 1single;
% udata = uout;

% double pulse
load 1double1a;
udata = ud_out;

% interpolate onto a larger grid
N = 512;
L = -xout(1);
L = 200;
[xout, uout] = fsolveequation(x, udata, par, N, L, config);
udata = uout;

% extract speed and wave
par.c = udata(end);
uwave = udata(1:end-1);

% differentiation matrices
Fourier = strcmp(config.method,'Fourier');
if Fourier
    [D, D2, D3, D4, D5] = D_fourier(N, L);
else
    % grid spacing
    h = xout(2) - xout(1);
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);
end

% check to see that our pulse (numerically) solves the eq
% and that derivative is eigenvector of J with eigenvalue 0
[F,J] = integratedequation(uwave,par,N,config,D,D2,D3,D4,D5);
% plot(xout, F)
% plot(xout, J*D*uout)

% % eigenvalues of integrated equation
% num = 50;
% center = -100;
% [int_lambda, ~, ~] = eigs_linear(xout, uwave, par, config, num, center, 'integrated');

% eigenvalues of nonintegrated equation
exp_wt = 0.001;
% [lambda_eig, V_eig, ~] = eig_linear(xout, uwave, par, config, 'nonintegrated', exp_wt);
num    = 30;
center = 0.01;
[lambda, V, J] = eigs_linear(xout, uwave, par, config, num, center, 'nonintegrated', exp_wt);

figure;
plot(lambda, '.');
% axis([-1e-3 1e-3 -1e5 1e5]);

