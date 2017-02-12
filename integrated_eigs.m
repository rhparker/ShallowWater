%% eigenvalues of integrated operator

% % single pulse
% load uc2;
% index = 10;
% 
% udata = uc(1:end-1, index);
% par.c = uc(end,index);
% 
% xout  = x;
% uout  = udata;

% % interpolate onto a larger grid
% N = 2*N;
% [xout, uout] = shallow_fourier_interp(x, uc(:,index), b, L, N);
% udata = uout(1:end-1);

% load 1usingle;
% udata = uout;

load 1udouble2a;
udata = ud_out(1:end-1);

b = par.b;

% differentiation matrices
[D, D2, D3, D4, D5] = D_fourier(N, L);

% check to see that our pulse (numerically) solves the eq
% and that derivative is eigenvector of J with eigenvalue 0
[F,J] = integratedshallow(udata,b,par,N,D,D2,D3,D4);
% plot(xout, F)
% plot(xout, J*D*uout)

num = 50;
center = -1000;
[int_lambda, ~, ~] = eigs_linear(xout, udata, b, par, L, N, num, center, 'integrated');
% [lambda, V, J] = eig_linear(xout, udata, b, par, L, N, 'integrated');
[lambda, V, J] = eig_linear(xout, udata, b, par, L, N, 'nonintegrated');

% figure;
% plot(lambda, zeros(length(lambda)),'.');
% title('eigenvalues of E"(u*), single pulse') 
% axis([-20 20 -1 1]);

