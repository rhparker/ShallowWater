% eigenvalues of linearization about stationary wave
% uses eigs instead of eig for speed
function [lambda, V, J] = eigs_linear(x, u, par, config, num, center, version, exp_wt)
    % exponential weight, if given
    if ~exist('exp_wt','var')
        exp_wt = 0;
    end

    % get the Jacobian of our eq around u
    J = get_jacobian(x, u, par, config, version, exp_wt);
    
    % run eigs on J
    opts.tol = 10^(-10);
    [V, lambda_D] = eigs(J, num, center, opts);
    lambda = diag(lambda_D);
end
