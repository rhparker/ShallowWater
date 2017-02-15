%% integrated shallow water wave equation
% parameters are b and c

function [F,J] = integratedshallow(u,par,N,D,D2,D3,D4,D5)
% returns the right-hand side of our equation

% operator

% linear part
LN = (2/15)*D4 - par.b*D2 + sparse(1:N,[1:N],par.c,N,N);
% we will need first two derivatives of u, compute them here
Du  = D*u;
D2u = D2*u;
% operator includes linear and nonlinear terms
F = LN*u + (3/2)*u.*u + (1/2)*Du.*Du + u.*D2u;

% Jacobian
if nargout > 1 
    ConstTerms = sparse(1:N,[1:N],3*u,N,N) + sparse(1:N,[1:N],D2u,N,N);
    D1Terms    = sparse(1:N,[1:N],Du,N,N)*D;
    D2Terms    = sparse(1:N,[1:N],u,N,N) *D2;
    J = LN + ConstTerms + D1Terms + D2Terms;
end

end
