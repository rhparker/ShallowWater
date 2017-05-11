% load 0single;
% 
% N = length(xout);
% L = -xout(1)*2;

N = 256;
L = 25;

% generate chebyshev (non-uniform) grid
xc = L*chebdif(N,1);

% % since no longer periodic BCs, ditch the initial point
% xout = xout(2:end);
% uwave = uout(2:end-1);
% uwave = (uwave + flip(uwave))/2;
% par.c = ud_out(end);
% 
% % cubic spline interpolation onto the chebyshev grid
% u = spline(xout, uwave, xc);
par.c = 36/169;
u = (105/338)*sech( xc / (2 * sqrt(13) ) ).^4;

config.method = 'Chebyshev';
config.BC = 'Neumann';
uold = [u; par.c];

[D, D2, D3, D4, D5] = D_cheb(N, L, config);
[F, J] = integratedequation(xc,u,par,N,config,D,D2,D3,D4,D5);

iter = 1000;
[xnew, unew] = fsolveequation(xc, uold, par, N, L, config, iter);
