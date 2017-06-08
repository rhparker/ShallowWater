% convert Fourier spectral solution to Chebyshev solution

function [x_cheb, u_cheb] = fourtocheb(x_four, u_four, config)
    L = ceil(abs(x_four(1)));         % domain size
    N = length(x_four);               % number of grid points
    par.c = u_four(end);
    
    % change config to Chebyshev / Neumann BCs
    config.method = 'Chebyshev';
    config.BC = 'Neumann';
    
    % make domain not periodic by reflecting initial point to end
    N = N + 1;
    uin = [u_four(1:end-1) ; u_four(1)];
    xin = [x_four ; -x_four(1)];
    [~, ~, ~, ~, ~, x_cheb] = D_cheb(N, L, config);
    uin = [spline(xin, uin, x_cheb); par.c];
    
    [~, u_cheb] = fsolveequation(x_cheb, uin, par, N, L, config, 1000);
end