% 5th order KdV equation

% a is optional exponential weight for Jacobian

function [F,J] = KdV(x,u,par,config,D,D2,D3,D4,D5,a)
% returns the right-hand side of our equation

%% operator

N = length(x);
L = ceil(abs(x(1)));
h = 2*L/N;

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
    
    % symmetry parameter
    if isfield(config, 'symmetry')
        % enforce symmetry by requiring (discrete) integral of 
        % ( u(x) - u(-x) )^2 = 0, i.e. discrete L2 norm squared
        if strcmp(config.symmetry, 'L2squaredflip')
            % how much to weight the L2 difference integral
            weight = 1/max(F);

            if strcmp(config.BC, 'periodic')
                usum = u(2:N);
                symm_J   = weight * 4*h*[0 ; usum - flip(usum) ]; 
            else
                usum = u(1:N);
                symm_J   = weight * 4*h*[ usum - flip(usum) ]; 
            end

            % discrete integrals and Jacobians of them
            symm_sum = weight * h*sum( (usum - flip(usum)).^2 );
            F = [F; symm_sum];
            J = [J; symm_J'];
        end

    end

end
