% continuation code for domain size
% actually does not use continuation code...

load 100F_50;

% which double pulse to pick
uout = ud_out_2_512;
xout = xout_512;

% parameters
L = ceil(abs(xout(1)));         % domain size
N = length(xout);               % number of grid points
h = (2*L)/N;
par.c = uout(end);

delta = 10;
steps = 10;
output_x = [ xout ];
output_u = [ uout ];
output_L = [ L ];

for index = 1:steps
    L = L+delta;
    [x, u] = double_pulse(L);
    output_x = [ output_x x ];
    output_u = [ output_u u ];
    output_L = [ output_L L ];
end

save 100F_domain output_x output_L output_u config

