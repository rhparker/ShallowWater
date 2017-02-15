%% generate double pulses

% load data for single pulse
% load uc4;

load KdV_uc_Fourier;
% load KdV_uc_fdiff;

% which equation to use
shallow = strcmp(config.equation,'shallow');

index = 200;    % KdV Fourier c = 6.2848
index = 186;    % KdV fdiff Neumann

% wave data and speed c
uout  = uc(:, index);
uwave = uout(1:end-1);
par.c = uc(end,index);

xout  = x;

h = xout(2) - xout(1);        % grid spacing
N = length(xout);             % number of grid points

% % interpolate onto a larger grid, if desired
% N = 2*N;
% [xout, uout] = shallow_fourier_interp(x, uc(:,index), b, L, N);

%% make half-wave from full wave

xhalf = xout( N/2+1 : end );
uhalf = uwave( N/2+1 : end );

% find the spatial eigenvalues
if shallow
    % shallow water eq, roots of (2/15) nu^4 - b nu^2 + c == 0
    nu = roots([(2/15) 0 -par.b 0 par.c]);
else
    nu = roots([1 0 -1 0 par.c]);
end

% oscillations frequency is imag(nu)
% oscillations decay with constant re(nu)
decay = abs(real(nu(1)));
freq  = abs(imag(nu(1)));
% spacing between pulses is expected to be an integer multiple of this
spacing = pi/freq;

% % for check, make plot of oscillations at foot of pulse
% % and find their expected frequency and exp decay rate
% cutoff = 5;
% osc_plot(xhalf, uhalf, b, par.c, L, cutoff);

% % do the join based on the frequency we computed
% % find the first minimum, which is the bottom dip of the wave
% [pks, locs] = findpeaks(-uhalf);
% start_pt = locs(1);
% start_x  = xhalf(start_pt);

% % x value and point to do the join
% join_x  = start_x + (minmax - 1)*spacing;

% do the join by finding the min/max of half-wave directly
% scale solution exponentially by appropriate scale factor
% since easier to find min/max of that
uscaled = uhalf.*exp(decay*xhalf);
% use finite diff/Neumann BCs to find derivative of uscaled
Nhalf   = length(xhalf);
FD      = D_fdiff_Neumann(Nhalf, h);
Duhalf  = FD*uscaled./uscaled - decay;
% interpolate derivative on finer grid to find zeros
Nfine       = 10000;
xfine       = linspace(0,L,Nfine)';
Duhalf_fine = interp1(xhalf,Duhalf,xfine);

% zero derivative is where we have a sign change
zDer = find(diff(sign(Duhalf_fine)));
% but derivative is discontinuous, so we only only want every other one
numMinMax = 5;                      % how many min/max we want
zDer      = zDer(2*(1:numMinMax));  % take every other one, starting at second 
zDer_x    = xfine(zDer);            % x values of deriative

% which min/max we want
minmax = 1;

% x value for join
join_x = zDer_x(minmax);

% % add this line to find half-way waves
% join_x = (zDer_x(minmax) + zDer_x(minmax+1) )/ 2;

% % add this line with minmax = 1 to find pulse 0 (half way to first min)
% % we might not be able to do this for shallow water eq
% join_x = join_x / 2;

% where to join the waves
join_pt = round(join_x / h) + 1;

% right half-wave
center_pt  = N/2 + 1; 
ud_half    = uout( join_pt : center_pt + join_pt - 1 );
ud_right   = flip( ud_half(1:end-1) );
% for periodic BCs, remove the final point
ud_right   = ud_right(1:end-1);
ud         = [ ud_half ; ud_right; par.c ];

% run joined pulse through Newton solver
% Newton solver on right half wave
iter = 1000;
[~, ud_out] = fsolveequation(xout, ud, par, N, L, config, iter);

% plot double wave before and after Newton solver
figure;
% plot the half wave
plot(xout, ud(1:end-1), xout, ud_out(1:end-1));
legend('initial guess','Newton solver output');
title('double pulse');


