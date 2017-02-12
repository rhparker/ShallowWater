% eigenvalues of linearization about stationary wave
function [lambda, V, J] = eig_linear(x, u, b, par, L, N, version)
	% whether to use integrated version or not
    integrated = strcmp(version, 'integrated');
    
    % generate differentiation matrices
    [D, D2, D3, D4, D5] = D_fourier(N, L);
    if integrated
        [F,J] = integratedshallow(u,b,par,N,D,D2,D3,D4);
    else
        [F,J] = shallow(u,b,par,N,D,D2,D3,D4,D5); 
    end
    [V, LD] = eig(full(J));
    lambda = diag(LD);
end
