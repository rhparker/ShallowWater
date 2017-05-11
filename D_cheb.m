%% Chebyshev differentiation matrices
% for interval [-L, L]

function [D, D2, D3, D4, D5, x] = D_cheb(N, L, config)

% generate Chebyshev matrices
[x, DM] = chebdif(N, 5);
D  = DM(:,:,1);
D2 = DM(:,:,2);
D3 = DM(:,:,3);
D4 = DM(:,:,4);
D5 = DM(:,:,5);

% scale to interval [-L, L]
x = L*x;
scale = 1/L;
D  = D  * scale;
D2 = D2 * scale^2;
D3 = D3 * scale^3;
D4 = D4 * scale^4;
D5 = D5 * scale^5;

% Neumann BCs
if strcmp(config.BC, 'Neumann')
    D(1,:)  = zeros(1, N);
    D(N,:)  = zeros(1, N);
    D2(1,:) = zeros(1, N);
    D2(N,:) = zeros(1, N);
end

end