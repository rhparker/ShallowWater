%% eigenvalues of integrated operator

% % single pulse
% load 0single_fourier;
% load 0single_fdiff;
% load 1single;
% uout = u2048;
% xout = x2048;

% % double pulse
% load 2double1a;
% load 2double2a;
load 3double1a;
% load 4double1a;
uout = ud_out;

% default parameters
par.c = uout(end);              % wave speed
N = length(xout);               % current grid size
L = -xout(1);                   % current domain length

% % if we want, we can change paramaters

% % change grid size
% N = 10000;                  % finite diff
% N = 2048;                   % Fourier
% 
% % change domain size
% L = 100;
N = 1024;

if N ~= length(xout) || L ~= -xout(1)
    % interpolate onto a larger grid, or with a longer domain
    [xnew, unew] = fsolveequation(xout, uout, par, N, L, config, 1000);
else
    % if we don't interpolate
    xnew = xout; unew = uout;
end

% just the wave
uwave = unew(1:end-1);

% differentiation matrices
Fourier = strcmp(config.method,'Fourier');
if Fourier
    [D, D2, D3, D4, D5] = D_fourier(N, L);
else
    % grid spacing
    h = xnew(2) - xnew(1);
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);
end

% % if we want to use the zero solution, uncomment this
% uwave = zeros([length(uwave) 1]);

%
% several checks that we can do
%

% check to see that our pulse (numerically) solves the eq
% and that derivative is eigenvector of J with eigenvalue 0
% [F,J] = integratedequation(uwave,par,N,config,D,D2,D3,D4,D5);
% plot(xout, F)
% plot(xout, J*D*uout)

% % check eigenvalues of integrated equation
% num = 50;
% center = -100;
% [int_lambda, ~, ~] = eigs_linear(xout, uwave, par, config, num, center, 'integrated');
% plot(int_lambda, zeros(length(int_lambda)), '.');

%% eigenvalues

% eigenvalues of nonintegrated equation

% find the ideal exponential weight (at least I hope so)
% in terms of c
a = find_exp_wt(par.c);
a = .01;
% a = 0;

% eigenvalue of the constant solution for linearization about zero solution
lambda_const = -a^5 + a^3 - par.c * a;

% [lambda_eig, V_eig, ~] = eig_linear(xnew, uwave, par, config, 'nonintegrated', a);
num    = 50;
center = 0.1;
[lambda, V, J] = eigs_linear(xnew, uwave, par, config, num, center, 'nonintegrated', a);

figure;
plot(lambda, '.');
% title('weighted space eigs for single pulse, known solution, fdiff/Neumann')
% axis([-0.002 0.003 -0.05 0.05]);

% figure;
% plot(lambda, '.');
% title('weighted space eigs for single pulse, known solution, fdiff/Neumann')
% axis([-0.001 0.001 -0.01 0.01]);

% figure;
% index = 3;
% plot(xnew,V(:,index))
% title(strcat('eigenfunction corresponding to eigenvalue:  ',num2str(lambda(index))))

% % use for complex eigenfunctions
% indices = 4:5;
% lambdaV = lambda(indices);
% legendcell = cellstr(num2str(lambdaV))
% 
% figure;
% plot( xnew, real( V(:,indices)) );
% legend(legendcell);
% title('Real part of eigenfunctions')
% 
% figure;
% plot( xnew, imag( V(:,indices)) );
% legendcell = cellstr(num2str(lambdaV))
% legend(legendcell);
% title('Imaginary part of eigenfunctions')
