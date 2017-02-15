% 5th order KdV equation

function [F,J] = KdV(u,par,N,D,D2,D3,D4,D5)
% returns the right-hand side of our equation

%% operator

% 5th order KDV
LN = D5 - D3 + par.c*D;
F  = LN*u - 2*u.*(D*u);

%% Jacobian
if nargout > 1     
    % 5th order KDV
    J = LN - 2 * sparse(1:N,[1:N],D*u,N,N) - 2 * sparse(1:N,[1:N],u,N,N)*D;
end

end
