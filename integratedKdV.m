% integrated version of 5th order KdV equation

function [F,J] = integratedKdV(u,par,N,D,D2,D3,D4,D5)
% returns the right-hand side of our equation

%% operator

% 5th order KDV, integrated once
LN = D4 - D2 + sparse(1:N,[1:N],par.c,N,N);
F = LN*u - u.*u;

%% Jacobian
if nargout > 1 
    % 5th order KDV, intergrated once
    J = LN - 2 * sparse(1:N,[1:N],u,N,N);
end

end
