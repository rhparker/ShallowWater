%% shallow water wave equation
% parameters are b and c

function [F,J] = shallow(u,par,N,D,D2,D3,D4,D5,a)
% returns the right-hand side of our equation

% operator

% linear part
LN = (2/15)*D5 - par.b*D3 + sparse(1:N,[1:N],par.c,N,N)*D;
% we will need first three derivatives of u, compute them here
Du  = D*u;
D2u = D2*u;
D3u = D3*u;
% operator includes linear and nonlinear terms
F = LN*u + 3*u.*Du + 2*Du.*D2u + u.*D3u;

% Jacobian
if nargout > 1 
    if a~0
        [XD, XD2, XD3, XD4, XD5] = D_expwt(D, D2, D3, D4, D5, a);
        LN = (2/15)*XD5 - par.b*XD3 + sparse(1:N,[1:N],par.c,N,N)*XD;
        ConstTerms = sparse(1:N,[1:N],3*Du,N,N)   + sparse(1:N,[1:N],D3u,N,N);
        D1Terms    = sparse(1:N,[1:N],3*u,N,N)*XD  + sparse(1:N,[1:N],2*D2u,N,N)*XD; 
        D2Terms    = sparse(1:N,[1:N],2*Du,N,N) * XD2;
        D3Terms    = sparse(1:N,[1:N],u,N,N) * XD3;
        J = LN + ConstTerms + D1Terms + D2Terms + D3Terms;
    else
        ConstTerms = sparse(1:N,[1:N],3*Du,N,N)   + sparse(1:N,[1:N],D3u,N,N);
        D1Terms    = sparse(1:N,[1:N],3*u,N,N)*D  + sparse(1:N,[1:N],2*D2u,N,N)*D; 
        D2Terms    = sparse(1:N,[1:N],2*Du,N,N) * D2;
        D3Terms    = sparse(1:N,[1:N],u,N,N) * D3;
        J = LN + ConstTerms + D1Terms + D2Terms + D3Terms;
    end
end

end
