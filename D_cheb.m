%% Chebyshev differentiation matrices
% for interval [-L, L]

function [L1, L2, L3, L4, L5, x] = D_cheb(N, L, config)

% generate Chebyshev matrices
[x, DM] = chebdif(N, 5);

% for Dirichlet BCs, remove first and last point
% which is equiv to eliminating top/bottom row
% and right/left col
DM = DM(2:end-1,2:end-1,:);
x = L*x(2:end-1);

% individual diff matrices
D1 = DM(:,:,1);
D2 = DM(:,:,2);
D3 = DM(:,:,3);
D4 = DM(:,:,4);
D5 = DM(:,:,5);

% scale to interval [-L, L]
scale = 1/L;
D1 = D1  * scale;
D2 = D2 * scale^2;
D3 = D3 * scale^3;
D4 = D4 * scale^4;
D5 = D5 * scale^5;

% % for homogeneous Neumann BCs (on 1st derivative)
% % we take our interpolating polynomial to be of form
% % p(x) = (L-x)(L+x) q(x)
% 
% M1 = diag(L^2 - x.^2);
% M2 = diag( 1./ (L^2 - x.^2) );
% M3 = diag( x./ (L^2 - x.^2) );
% M4 = diag( x );
% 
% L1 = M1*D1*M2 -  2*M3;
% L2 = M1*D2*M2 -  4*M4*D1*M2 -  2*M2;
% L3 = M1*D3*M2 -  6*M4*D2*M2 -  6*D1*M2;
% L4 = M1*D4*M2 -  8*M4*D3*M2 - 12*D2*M2;
% L5 = M1*D5*M2 - 10*M4*D4*M2 - 20*D3*M2;

% integrated form (5th order)
% for homogeneous Neumann BCs (on 1st derivative)
% and right-sided Neumann BCs (on 2nd derivative)
% we take our interpolating polynomial to be of form
% p(x) = (L-x)^2 (L+x) q(x)

if isfield(config, 'form') && strcmp(config.form, 'nonintegrated')
    M0 = diag( L^3 - L^2*x - L*x.^2 + x.^3 );
    M1 = diag( -L^2 - 2*L*x + 3*x.^2);
    M2 = diag( -2*L + 6*x);
    M3 = 6*eye(length(x));
    Me = diag( 1 ./ (L^3 - L^2*x - L*x.^2 + x.^3) );
    
% nonintegrated form (4th order)
% for homogeneous Neumann BCs (on 1st derivative)
% we take our interpolating polynomial to be of form
% p(x) = (L-x)(L+x) q(x)
else
    M0 = diag(L^2 - x.^2);
    M1 = diag(-2*x);
    M2 = -2 * eye(length(x));
    M3 =  0 * eye(length(x));
    Me = diag( 1./ (L^2 - x.^2) );
end

L1 = M0*D1*Me +   M1*Me;
L2 = M0*D2*Me + 2*M1*D1*Me +    M2*Me;
L3 = M0*D3*Me + 3*M1*D2*Me +  3*M2*D1*Me +    M3*Me;
L4 = M0*D4*Me + 4*M1*D3*Me +  6*M2*D2*Me +  4*M3*D1*Me;
L5 = M0*D5*Me + 5*M1*D4*Me + 10*M2*D3*Me + 10*M3*D2*Me;

end