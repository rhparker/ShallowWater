function [F,J] = integratedequation(x,u,par,N,config,D,D2,D3,D4,D5,usymm)

if ~exist('usymm','var')
    usymm = u;
end

% if specified, use shallow water equation
if strcmp(config.equation,'shallow')
    [F,J] = integratedshallow(u,par,N,D,D2,D3,D4,D5);
    
% otherwise use 5th order KdV equation
else
    [F,J] = integratedKdV(x,u,par,config,D,D2,D3,D4,D5,usymm);
end
