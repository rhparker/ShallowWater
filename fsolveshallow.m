function [xout, uout] = fsolveshallow(xold, uold, b, N, L, iter, c)
%% setup

% max number of iterations to use (default is 100)
if ~exist('iter','var')
    iter = 100;
end

% if we don't specify a new c, then take it from the 
% last element of uin
if ~exist('c','var')
    par.c = uold(end);
else
    par.c = c;
end

% extract wave data
u = uold(1:end-1);

% if N is different from size of grid xold, then
% interpolate onto a new grid
if N ~= length(xold)
    % if we have a nonperiodic domain, i.e. [-L, L]
    if xold(1) + xold(end) == 0
        xout = linspace(xold(1), xold(end),N)';
    else
        xout = linspace(xold(1),-xold(1),N+1)';
        xout = xout(1:end-1);
    end
    u = interp1(xold,uold(1:end-1),xout);

    % since we have no final grid point (periodic)
    % will have some NaN values at the end of u
    % replace these with 0
    u(isnan(u)) = 0;
else
    xout = xold;
end

% differentiation matrices
[D, D2, D3, D4, D5] = D_fourier(N, L);


%% solve nonlinear problem using fsolve

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',iter);

% call fsolve
[uout,fval] = fsolve(@(u) integratedshallow(u,b,par,N,D,D2,D3,D4),u,options);

% reappend c to the output vector
uout = [uout ; par.c];
