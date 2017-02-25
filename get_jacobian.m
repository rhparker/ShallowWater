% get Jacobian for equation (integrated or regular) about solution u
% optional parameter exp_wt for exponential weighted space

function J = get_jacobian(x, u, par, config, version, exp_wt)
    % number of grid points is length of x
    N = length(x);
    
    % generate differentiation matrices
    Fourier = strcmp(config.method,'Fourier');
    if Fourier
        % length of domain is negative of 1st element of x
        % since we are always on [-L, L]
        L = -x(1);
        [D, D2, D3, D4, D5] = D_fourier(N, L);
    else
        % grid spacing: dist between 1st two elements of x
        h = x(2) - x(1);
        [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);
    end

    % if no exponential weight supplied, it is zero (unweighted)
    if ~exist('exp_wt','var')
        exp_wt = 0;
    end

    % whether to use integrated version or not
    integrated = strcmp(version, 'integrated');
    if integrated
        [~,J] = integratedequation(u,par,N,config,D,D2,D3,D4,D5);
    else
        [~,J] = equation(u,par,N,config,D,D2,D3,D4,D5,exp_wt); 
    end
end