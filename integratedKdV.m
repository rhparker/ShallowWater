% integrated version of 5th order KdV equation

function [F,J] = integratedKdV(x,u,par,config,D,D2,D3,D4,D5,usymm)
% returns the right-hand side of our equation

%% operator

N = length(x);
L = ceil(abs(x(1)));
h = 2*L/N;

% 5th order KDV, integrated once
LN = D4 - D2 + sparse(1:N,[1:N],par.c,N,N);
F = LN*u - u.*u;

%% Jacobian
if nargout > 1 
    % 5th order KDV, intergrated once
    J = LN - 2 * sparse(1:N,[1:N],u,N,N);
end

% symmetry parameter
if isfield(config, 'symmetry')

    % enforce symmetry by comparing to previous step
    % assumption is that usymm is symmetric (even)
    if strcmp(config.symmetry, 'compare')
        if strcmp(config.BC, 'periodic')
            lowerb = 1;
            upperb = N;
        end

        % derivative and 2nd derivative of u
        Du  = D*u;
        D2u = D2*u;
        
        weight = 1;
        % discrete integrals and Jacobians of them
        symm_sum = weight * h * sum( Du.*(u - usymm) );
        symm_J   = weight * h * ( D2u.*Du + Du.*Du - D2u.*usymm );
        F = [F; symm_sum];
        J = [J; symm_J'];

    % enforce symmetry by requiring (discrete) integral of 
    % ( u(x) - u(-x) )^2 = 0, i.e. discrete L2 norm squared
    elseif strcmp(config.symmetry, 'L2squaredflip')
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

    % enforce symmetry by requiring discrete L2 norm of  
    % u(x) - u(-x) be zero
    elseif strcmp(config.symmetry, 'L2flip')
        if strcmp(config.BC, 'periodic')
            usum = u(2:N);
        end
        % discrete integrals and Jacobians of them
        symm_sum_sq = sum( (usum - flip(usum)).^2 );
        symm_sum    = sqrt(h) * sqrt( symm_sum_sq );
        symm_J      = (2 * sqrt(h) / symm_sum) * [0 ; usum - flip(usum) ]; 
        F = [F; symm_sum];
        J = [J; symm_J'];
    end
end

end
