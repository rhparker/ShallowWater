% continuation code for domain size
% actually does not use continuation code...
% parameters

load 1500F_10;

xout = xout_1024;

L = ceil(abs(xout(1)));         % domain size
N = length(xout);               % number of grid points
h = (2*L)/N;
par.c = uout(end);

delta = 5;
steps = 5;
output_x = [];
output_u = [];
output_L = [];

L = 105;

for index = 1:steps
    [x, u] = double_pulse(L);
    output_x = [ output_x x ];
    output_u = [ output_u u ];
    output_L = [ output_L L ];
    L = L+delta;
end

save 1500F_domain_105 output_x output_L output_u config

