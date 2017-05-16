function [F,J] = equation(x,u,par,N,config,D,D2,D3,D4,D5,a)

% if exponential weight is not specified, pass on 0
if ~exist('a','var')
    a = 0;
end

% if specified, use shallow water equation
if strcmp(config.equation,'shallow')
    [F,J] = shallow(x,u,par,N,D,D2,D3,D4,D5,a);
    
% otherwise use 5th order KdV equation
else
    [F,J] = KdV(x,u,par,config,D,D2,D3,D4,D5,a);
end 
    