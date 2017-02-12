%% Fourier differentiation matrices
% for interval [-L, L]

function [D, D2, D3, D4, D5] = D_fourier(N,L)

% generate Fourier matrices
[~,D]  = fourdif(N,1);
[~,D2] = fourdif(N,2);
[~,D3] = fourdif(N,3);
[~,D4] = fourdif(N,4);
[~,D5] = fourdif(N,5);

% scale to interval [-L, L]
scale = pi/L;
D  = D  * scale;
D2 = D2 * scale^2;
D3 = D3 * scale^3;
D4 = D4 * scale^4;
D5 = D5 * scale^5;

end