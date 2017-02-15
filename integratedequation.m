function [F,J] = integratedequation(u,par,N,config,D,D2,D3,D4,D5)

% if specified, use shallow water equation
if strcmp(config.equation,'shallow')
    [F,J] = integratedshallow(u,par,N,D,D2,D3,D4,D5);
    
% otherwise use 5th order KdV equation
else
    [F,J] = integratedKdV(u,par,N,D,D2,D3,D4,D5);
    
end
    
