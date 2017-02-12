% eigenvalues of linearization about stationary wave
% uses eigs instead of eig for speed
function [lambda, V, J] = eigs_linear(x, u, b, par, L, N, num, center, version)
    % whether to use integrated version or not
    integrated = strcmp(version, 'integrated');

    % generate differentiation matrices
    [D, D2, D3, D4, D5] = D_fourier(N, L);
    if integrated
        [F,J] = integratedshallow(u,b,par,N,D,D2,D3,D4);
    else
        [F,J] = shallow(u,b,par,N,D,D2,D3,D4,D5); 
    end
    opts.tol = 10^(-10);
    [V, lambda_D] = eigs(J, num, center, opts);
    lambda = diag(lambda_D);
end
