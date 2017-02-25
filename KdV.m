% 5th order KdV equation

% a is optional exponential weight for Jacobian

function [F,J] = KdV(u,par,N,D,D2,D3,D4,D5,a)
% returns the right-hand side of our equation

%% operator

% 5th order KDV
LN = D5 - D3 + par.c*D;
F  = LN*u - 2*u.*(D*u);

%% Jacobian
if nargout > 1     
    % 5th order KDV
    % exponentially weighted version
    if a ~= 0
        [XD, XD2, XD3, XD4, XD5] = D_expwt(D, D2, D3, D4, D5, a);
        LN = XD5 - XD3 + par.c*XD;
        J = LN - 2 * sparse(1:N,[1:N],D*u,N,N) - 2 * sparse(1:N,[1:N],u,N,N)*XD;
    else
    % normal version
        J = LN - 2 * sparse(1:N,[1:N],D*u,N,N) - 2 * sparse(1:N,[1:N],u,N,N)*D;
    end

end
