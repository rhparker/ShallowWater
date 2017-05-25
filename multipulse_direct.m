%% generate multiple pulses

% load from continuation data

load ucKdV_Fourier_256;
% index = 234;        % closest to 2.5
% index = 451;        % closest to 5
% index = 652;        % closest to 7.5
index = 844;        % closest to 10

% load ucKdV_Cheb_256;
% index = 189;        % closest to 2.5
% index = 365;        % closest to 5
% index = 528;        % closest to 7.5
% index = 684;        % closest to 10

% wave data and speed c
uout  = uc(:, index);
xout  = x;

% load pregenerated single pulse

% load 100C;

% extract speed parameter
par.c = uout(end);

% which equation to use
shallow = strcmp(config.equation,'shallow');

% extract length and grid point parameters
L = ceil(abs(xout(1)));         % domain size
N = length(xout);               % number of grid points
h = (2*L)/N;
uwave = uout(1:end-1);

% since using Dirichlet BCs on Chebyshev
if (strcmp(config.method, 'Chebyshev'))
    N = N + 2;
end

% adjust c if we want (to standardize)
 
par.c = 10;
uout(end) = par.c;

% adjust N or L if we want to

% N = 1024;
% L = 50;
% h=2*L/N;

% symmetry_after = false;
symmetry_after = true;

if symmetry_after
    config.symmetry = 'none';
else
    config.symmetry = 'L2squaredflip';
end

% if we change anything, need to send through fsolve again
[xout, uout] = fsolveequation(x, uout, par, N, L, config);
uwave = uout(1:end-1);
h = (2*L)/N;

% chebychev points listed backwards, so flip around 
% for double pulse construction
if strcmp(config.method,'Chebyshev')
    xout = flip(xout);
    uout = [ flip(uout(1:end-1)) ; par.c ];
end

findsymm(xout, uout, config);

%% make half-wave from full wave

% where to split for the half wave depends on whether
% we have periodic BCs or not
if strcmp(config.BC,'periodic')
    center = N/2 + 1;
else
    if mod(N,2) == 0
        center = N/2 + 1;
    else
        center = (N+1)/2;
    end
end

% if we use Chebyshev, interpolate onto a uniform grid
if strcmp(config.method,'Chebyshev')
    x_uniform = linspace(-L, L, N)';
    xhalf = x_uniform( center : end );
    uout = [interp1(xout, uwave, x_uniform) ; par.c];
else
    xhalf = xout( center : end );
end

% half wave
uwave = uout(1:end-1);
uhalf = uwave( center : end );

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

% linear interpolation
Duhalf_fine = interp1(xhalf,Duhalf,xfine);
Duhalf_fine(isnan(Duhalf_fine)) = 0;

% zero derivative is where we have a sign change
zDer = find(diff(sign(Duhalf_fine)));
% but derivative is discontinuous, so we only only want every other one
numMinMax = 2;                      % how many min/max we want
zDer      = zDer(2*(1:numMinMax));  % take every other one, starting at second 
zDer_x    = xfine(zDer);            % x values of deriative

% find join points using known spacing
start = 1;
index = 2;
join_x = zDer_x(start) + (index - 1)*(spacing/2);

% where to join the waves
join_pt = round(join_x / h) - 1;

% half waves and middle wave
u_left    = uout( 1 : center + join_pt - 1 );
u_right   = flip( u_left(1:end-1) );
u_middle  = uout( center - join_pt + 2: center + join_pt - 1);

num_waves = 2;
um_center = repmat(u_middle, num_waves-1, 1);

% construct multiple pulses
um = [ u_left ; um_center; u_right ];

% output too long, so need to chop off the ends
N_m = length(um);
if mod(N_m, 2) == 0
    center_m = N_m / 2 + 1;
else
    center_m = (N_m + 1)/2;
end

center_diff = center_m - center;
if mod(N_out, 2) == 0
    um = um( center_diff + 1 : center_diff + N );
else
    um = um( center_diff + 1 : center_diff + N );
end
% if we use Chebyshev, interpolate back onto Chebyshev grid
if strcmp(config.method,'Chebyshev')
    um = interp1(x_uniform, um, xout);
 
    % flip back around 
    xout = flip(xout);  
    um = flip(um);
end

% append speed c to end

um = [um ; par.c];

% run joined pulse through Newton solver
% Newton solver on right half wave
iter = 1000;
[~, um_out] = fsolveequation(xout, um, par, N, L, config, iter);

% if we enforce symmetry afterwards, send back through Newton
% solver with symmetry enforced

if symmetry_after
    % save this in case we want to compare
    um_old = um_out;
    % enforce symmetry
    config.symmetry = 'L2squaredflip';
    [~, um_out] = fsolveequation(xout, um_out, par, N, L, config, iter);
end

% plot double wave before and after Newton solver
figure;
plot(xout, um(1:end-1), xout, um_out(1:end-1));

legend('initial guess','Newton solver output');
title('double pulse');

findsymm(xout, um_out, config);


