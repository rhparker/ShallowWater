% eigenvalues of linearization about stationary wave
function [lambda, V, J] = eig_linear(x, u, par, config, version, exp_wt)
    % exponential weight, if given
    if ~exist('exp_wt','var')
        exp_wt = 0;
    end

    % get the Jacobian of our eq around u
    J = get_jacobian(x, u, par, config, version, exp_wt);
    
    % run eig on this; requires a full matrix
    [V, LD] = eig(full(J));
    lambda = diag(LD);
end
