function [F,J] = equation(u,par,N,config,D,D2,D3,D4,D5,a)

% if specified, use shallow water equation
if strcmp(config.equation,'shallow')
    [F,J] = shallow(u,par,N,D,D2,D3,D4,D5);
    
% otherwise use 5th order KdV equation
else
    [F,J] = KdV(u,par,N,D,D2,D3,D4,D5,a);
    
end 
    