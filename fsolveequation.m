function [xout, uout] = fsolveequation(xold, uold, par, N, L, config, iter)
%% setup

% max number of iterations to use (default is 100)
if ~exist('iter','var')
    iter = 100;
end

% extract wave data, N, and L 
u     = uold(1:end-1);
N_old = length(xold);
L_old = abs( xold(1) );

% differentiation matrices
if strcmp(config.method,'Fourier')
    [D, D2, D3, D4, D5] = D_fourier(N, L);
elseif strcmp(config.method,'Chebyshev')
    [D, D2, D3, D4, D5, xout] = D_cheb(N+2, L, config);
% finite differences
else
    % grid spacing
    h = xout(2) - xout(1);
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);
end

% if N is different from size of grid xold, 
% or L is different from domain of xold,
% then interpolate onto a new grid
if (N ~= N_old) || (L ~= L_old)
    % if we have a nonperiodic domain, i.e. [-L, L]
    if xold(1) + xold(end) == 0
        if strcmp(config.method,'Chebyshev')
            % compute xout when we compute diff matrices
        else
            % xout = linspace(xold(1), xold(end),N)';
            xout = linspace(-L, L, N)';
        end
    % for periodic domain, have to remove final point
    else
        % xout = linspace(xold(1),-xold(1),N+1)';
        xout = linspace(-L, L, N+1)';
        xout = xout(1:end-1);
    end
    
    if strcmp(config.method,'Fourier')
        % fourier interpolation
        u = interpft(uold(1:end-1), N);
    else
        % linear interpolation
        u = interp1(xold,uold(1:end-1),xout);
        % will have some NaN values at the ends of u
        % replace these with 0
        u(isnan(u)) = 0;
    end

else
    xout = xold;
end

%% solve nonlinear problem using fsolve

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',iter);
options.TolFun = 1e-16;
options.TolX = 1e-16;
% call fsolve

if isfield(config, 'form') && strcmp(config.form, 'nonintegrated');
    [uout,fval] = fsolve(@(u) equation(xout,u,par,N,config,D,D2,D3,D4,D5),u,options);
else
    [uout,fval] = fsolve(@(u) integratedequation(xout,u,par,N,config,D,D2,D3,D4,D5),u,options);
end
% [uout,fval] = fsolve(@(u) equation(u,par,N,config,D,D2,D3,D4,D5),u,options);

% reappend c to the output vector
uout = [uout ; par.c];
