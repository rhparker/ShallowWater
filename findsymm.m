% find symmetry of function
function [xout, uout] = findsymm(xin, uin, config)

% grid points to use for interpolation
N = 10000;

L = -xin(1);
xout = linspace(-L, L, N+1)';
xout = xout(1:end-1);
uout = interp1(xin,uin,xout);
uout(isnan(uout)) = 0;


% % for periodic BCs, remove the first point
% if strcmp(config.BC, 'periodic')
%     xout = xout(2:end);
%     uout  = uout(2:end);
% end


end