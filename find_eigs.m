%% eigenvalues of integrated operator

% single pulse
% load 0single_fourier;
% load 0single_fdiff;
% load 5single;

load 5double1;
uout = ud_out;

% default parameters
par.c = uout(end);              % wave speed
N = length(xout);               % current grid size
L = -xout(1);                   % current domain length

N=1024;

% % if we want, we can change paramaters

% % change grid size
% N = 20000;                  % finite diff
% N = 2048;                   % Fourier
% N = 128;
% N = 4096;
% N = 400;
N = 10000;

% % change domain size
% L = 100;
% L = 25;
% L = 50;

% % change speed c (to standardize between methods)
% par.c = 82.5;

% if we change stuff, run through Newton solver
if N ~= length(xout) || L ~= -xout(1) || par.c ~= uout(end)
    % interpolate onto a larger grid, or with a longer domain, or with
    % different c
    [xnew, unew] = fsolveequation(xout, uout, par, N, L, config, 10000);
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
    h = 2*L/N;
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC, 3);
end

% % if we want to use the zero solution, uncomment this
% uwave = zeros([length(uwave) 1]);

%
% several checks that we can do
%

% % check to see that our pulse (numerically) solves the eq
% % and that derivative is eigenvector of J with eigenvalue 0
[F,J] = equation(uwave,par,N,config,D,D2,D3,D4,D5);
% [F,J] = integratedequation(uwave,par,N,config,D,D2,D3,D4,D5);
% plot(xnew, F)
% plot(xnew, J*D*uwave)

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
% a = a / 2;
a = 0
% a = .2;

% eigenvalue of the constant solution for linearization about zero solution
lambda_const = -a^5 + a^3 - par.c * a;

num    = 100;
% center = 1.75i;
% center = -0.004+0.4i;
center = 3;

% % use eig
% [lambda, V, J] = eig_linear(xnew, uwave, par, config, 'nonintegrated', a);
% % use eig
[lambda, V, J] = eigs_linear(xnew, uwave, par, config, num, center, 'nonintegrated', a);

eig_plot = true;

if eig_plot
    figure;
    plot(lambda, '.');
    plot_title   = ['Eigenvalues using eigs of double pulse 1, exp wt a = ',num2str(a)];
    method_title = [config.method,', ',config.BC,' (N = ',num2str(N),', L = ',num2str(L),') '];
    % eig_title  = ['lambda = ',num2str(lambda(1))];
    title({plot_title, [method_title, ' ']}); 
end

% for exponentially weighted space, if we have a good
% separation, we can extract the eigenvalues which are
% to the R of the essential spectrum
if a ~ 0
    cutoff  = -2;
    indices = find(real(lambda) > cutoff);
    eVals = lambda(indices);
    eVecs = V(:, indices);
    
    % for eigenvalues close to imag axis, we can fsolve them to there
    index = 1;
    [vout, lout] = eig_solve(J, i*imag(eVals(index)), eVecs(:,index), 'fix');
    max_before = max( abs( J*eVecs(:,index) - eVals(index)*eVecs(:,index)) );
    max_after  = max( abs( J*vout - lout*vout));
    
% otherwise grab the eigenvalues off the real axis
else
    cutoff = 0.1;
    indices = find(abs(real(lambda)) > cutoff);
    eVals = lambda(indices);
    eVecs = V(:, indices);
end

% plot(xnew, exp(-a*xnew).*real(vout));
% title('real part of eigenfunction, eigenvalue 0.6423i');

% plot(xnew, exp(-a*xnew).*real(vout));
% title('real part of unweighted eigenfunction, eigenvalue 0.6423i');

% plot(xnew, imag(vout));
% title('imag part of eigenfunction, eigenvalue 0.6423i');


% figure;
% plot(lambda, '.');
% plot_title   = ['Zoom of eigenvalues using eigs of double pulse 2, exp wt a = ',num2str(a)];
% method_title = [config.method,', ',config.BC,' (N = ',num2str(N),', L = ',num2str(L),') '];
% eig_title  = ['lambda = ',num2str(lambda(1))];
% title({plot_title, [method_title, ' ']}); 
% axis([-0.4 0.4 -3 3]);

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
