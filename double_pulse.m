%% generate double pulses

% load from continuation data

load ucKdV_Fourier_256;
% load ucKdV_Fourier_128;

% which equation to use
shallow = strcmp(config.equation,'shallow');

% index = 200;    % KdV Fourier c = 6.2848;
% index = 186;    % KdV fdiff Neumann, c = 6.2755
% index = 500;    % KdV fdiff Neumann 2, c = 18.7068
% index = 266;    % KdV Fourier 2, c = 18.6714
% index = 232;    % KdV Fourier 2, c = 16.0357
% index = 436;    % KdVfdiff 2, c = 16.0390
% index = 214;
% index = 1000;     % KdVfdiff1000, c = 40.9355
% index = 537;      % KdVfourier128, c = 40.9572
% index = 591;      % KdVfourier256_40, c = 40.9273

index = 1465;       % KdV_Fourier_256, c = 40.9410
index = 1600;
index = 200;

% wave data and speed c
uout  = uc(:, index);
uwave = uout(1:end-1);
par.c = uc(end,index);
xout  = x;

L = -xout(1);                 % domain size
N = length(xout);             % number of grid points
h = (2*L)/N;
uwave = uout(1:end-1);

% adjust c if we want (to standardize)
% par.c = 6.275;
% par.c = 18.67;
% par.c = 16;
% par.c = 5;
% par.c = 40.9355;
% uout(end) = par.c;

% adjust N if we want to
% N = 1024;

% N=2000;
% L=50;

config.symmetry = 'L2squaredflip';

% if we change anything, need to send through fsolve again
[xout, uout] = fsolveequation(x, uout, par, N, L, config);
uwave = uout(1:end-1);
h = (2*L)/N;

findsymm(xout, uout, config)

%% make half-wave from full wave

% where to split for the half wave depends on whether
% we have periodic BCs or not
if strcmp(config.BC,'Neumann')
    offset = 0;
% periodic BCs
else
    offset = 1;
end

xhalf = xout(  N/2+offset : end );
uhalf = uwave( N/2+offset : end );

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

% do the join by finding the min/max of half-wave directly
% scale solution exponentially by appropriate scale factor
% since easier to find min/max of that
uscaled = uhalf.*exp(decay*xhalf);
% use finite diff/Neumann BCs to find derivative of uscaled
Nhalf   = length(xhalf);
FD      = D_fdiff_Neumann(Nhalf, h);
Duhalf  = FD*uscaled./uscaled - decay;
% interpolate derivative on finer grid to find zeros
Nfine       = 10001;
xfine       = linspace(0,L,Nfine)';
Duhalf_fine = interp1(xhalf,Duhalf,xfine);
Duhalf_fine(isnan(Duhalf_fine)) = 0;

% zero derivative is where we have a sign change
zDer = find(diff(sign(Duhalf_fine)));
% but derivative is discontinuous, so we only only want every other one
numMinMax = 3;                      % how many min/max we want
zDer      = zDer(2*(1:numMinMax));  % take every other one, starting at second 
zDer_x    = xfine(zDer);            % x values of deriative

% % which min/max we want
% minmax = 1;
% % x value for join
% join_x = zDer_x(minmax);

% another way to do this, using known spacing
start = 1;
index = 2;
join_x = zDer_x(start) + (index - 1)*(spacing/2);

% add this line to find half-way waves
% join_x = (zDer_x(minmax) + zDer_x(minmax+1) )/ 2;

% % add this line with minmax = 1 to find pulse 0 (half way to first min)
% % we might not be able to do this for shallow water eq
% join_x = join_x / 2;

% where to join the waves
join_pt = round(join_x / h)+1;

% right half-wave
center_pt  = N/2 + 1; 
ud_half    = uout( join_pt : center_pt + join_pt - 1 );
ud_right   = flip( ud_half(1:end-1) );
% for periodic BCs, remove the final point
ud_right   = ud_right(1:end-1);
ud         = [ ud_half ; ud_right; par.c ];

% run joined pulse through Newton solver
% Newton solver on right half wave
iter = 10000;
[~, ud_out] = fsolveequation(xout, ud, par, N, L, config, iter);

% plot double wave before and after Newton solver
figure;
% plot the half wave
plot(xout, ud(1:end-1), xout, ud_out(1:end-1));
legend('initial guess','Newton solver output');
title('double pulse');

findsymm(xout, ud_out, config);


