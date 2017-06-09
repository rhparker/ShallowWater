% continuation code for domain size
% actually does not use continuation code...
% parameters

load 100F;

L = ceil(abs(xout(1)));         % domain size
N = length(xout);               % number of grid points
h = (2*L)/N;
par.c = uout(end);

delta = 10;
steps = 5;
output_x = [];
output_u = [];
output_L = [];

L = 125;

for index = 1:steps
    [x, u] = double_pulse(L);
    output_x = [ output_x x ];
    output_u = [ output_u u ];
    output_L = [ output_L L ];
    L = L+delta;
end

save 100F_domain output_x output_L output_u config

