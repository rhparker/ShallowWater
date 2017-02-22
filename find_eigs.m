%% eigenvalues of integrated operator

% % single pulse
% load 1single;

% double pulse
% load 2double1a;
% load 2double2a;
load 1double1a;
uout = ud_out;

% wave speed
par.c = uout(end);

% interpolate onto a larger grid, or with a longer domain
N = length(xout);               % current grid size
% change grid size
% N = 10000;
N = 1024;

L = -xout(1);                   % current domain length
% change domain size
L = 100;

[xnew, unew] = fsolveequation(xout, uout, par, N, L, config);
uwave = unew(1:end-1);

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

%% eigenvalues

% eigenvalues of nonintegrated equation
exp_wt = 0.0001;
% [lambda_eig, V_eig, ~] = eig_linear(xout, uwave, par, config, 'nonintegrated', exp_wt);
num    = 50;
center = 0.001;
[lambda, V, J] = eigs_linear(xnew, uwave, par, config, num, center, 'nonintegrated', exp_wt);

figure;
plot(lambda, '.');
title('zoom of eigs centered near 0 for 2nd double pulse, N = 100000, L = 100')
axis([-0.002 0.003 -0.05 0.05]);

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
