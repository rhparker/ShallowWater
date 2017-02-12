%% shallow water wave equation
% parameters are b and c
% since we will keep b fixed and vary c by continuation, store c
% in struct as par.c

function [F,J] = shallow(u,b,par,N,D,D2,D3,D4,D5)
% returns the right-hand side of our equation

% operator

% linear part
LN = (2/15)*D5 - b*D3 + sparse(1:N,[1:N],par.c,N,N)*D;
% we will need first three derivatives of u, compute them here
Du  = D*u;
D2u = D2*u;
D3u = D3*u;
% operator includes linear and nonlinear terms
F = LN*u + 3*u.*Du + 2*Du.*D2u + u.*D3u;

% Jacobian
if nargout > 1 
    ConstTerms = sparse(1:N,[1:N],3*Du,N,N)   + sparse(1:N,[1:N],D3u,N,N);
    D1Terms    = sparse(1:N,[1:N],3*u,N,N)*D  + sparse(1:N,[1:N],2*D2u,N,N)*D; 
    D2Terms    = sparse(1:N,[1:N],2*Du,N,N) * D2;
    D3Terms    = sparse(1:N,[1:N],u,N,N) * D3;
    J = LN + ConstTerms + D1Terms + D2Terms + D3Terms;
end

end
