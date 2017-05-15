function [F,J] = integratedequation(x,u,par,N,config,D,D2,D3,D4,D5,usymm)

if ~exist('usymm','var')
    usymm = u;
end

% For Dirichlet BCs, we only pass interior points and 
% set ends to 0
if isfield(config, 'Dirichlet') && strcmp(config.Dirichlet, 'true')
    u = [0 ; u ; 0];
end

% if specified, use shallow water equation
if strcmp(config.equation,'shallow')
    [F,J] = integratedshallow(u,par,N,D,D2,D3,D4,D5);
    
% otherwise use 5th order KdV equation
else
    [F,J] = integratedKdV(x,u,par,config,D,D2,D3,D4,D5,usymm);
end

% For Dirichlet BCs, the endpoints are not variables, so
% the Jacobian does not depend on them, so we need to
% remove the first and last columns
if isfield(config, 'Dirichlet') && strcmp(config.Dirichlet, 'true')
    J = J(:, 2:end-1);
end

