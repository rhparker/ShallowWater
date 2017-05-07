function [xout, uout] = fsolveequation(xold, uold, par, N, L, config, iter)
%% setup

% max number of iterations to use (default is 100)
if ~exist('iter','var')
    iter = 100;
end

% which method to use
Fourier = strcmp(config.method,'Fourier');

% % if we don't specify a new c, then take it from the 
% % last element of uin
% if ~exist('c','var')
%     par.c = uold(end);
% else
%     par.c = c;
% end

% extract wave data, N, and L 
u     = uold(1:end-1);
N_old = length(xold);
L_old = -xold(1);

% if N is different from size of grid xold, 
% or L is different from domain of xold,
% then interpolate onto a new grid
if (N ~= N_old) || (L ~= L_old)
    % if we have a nonperiodic domain, i.e. [-L, L]
    if xold(1) + xold(end) == 0
        % xout = linspace(xold(1), xold(end),N)';
        xout = linspace(-L, L, N)';
    else
        % xout = linspace(xold(1),-xold(1),N+1)';
        xout = linspace(-L, L, N+1)';
        xout = xout(1:end-1);
    end
    
    % linear interpolation
%     u = interp1(xold,uold(1:end-1),xout);
%     % will have some NaN values at the ends of u
%     % replace these with 0
%     u(isnan(u)) = 0;

    % fourier interpolation
    u = interpft(uold(1:end-1), N);

else
    xout = xold;
end

% differentiation matrices
if Fourier
    [D, D2, D3, D4, D5] = D_fourier(N, L);
else
    % grid spacing
    h = xout(2) - xout(1);
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);
end


%% solve nonlinear problem using fsolve

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',iter);
options.TolFun = 1e-16;
options.TolX = 1e-16;
% call fsolve
[uout,fval] = fsolve(@(u) integratedequation(xout,u,par,N,config,D,D2,D3,D4,D5),u,options);
% [uout,fval] = fsolve(@(u) equation(u,par,N,config,D,D2,D3,D4,D5),u,options);

% reappend c to the output vector
uout = [uout ; par.c];
