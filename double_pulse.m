%% generate double pulses

function [x, u] = double_pulse(L_new)

% load from continuation data

% load ucKdV_Fourier_256;
% index = 234;        % closest to 2.5
% index = 451;        % closest to 5
% index = 652;        % closest to 7.5
% index = 844;        % closest to 10

% load ucKdV_Cheb_256;
% index = 189;        % closest to 2.5
% index = 365;        % closest to 5
% index = 528;        % closest to 7.5
% index = 684;        % closest to 10

% load ucKdV_Fourier_256_50;
% index = 502;          % closest to 10

% load ucKdV_Fourier_512_200;
% index = 456;          % closest to 10

% load 25C;
% load 100F;

load 100F_200;
xout = xout_1024;
uout = uout_1024;

% which equation to use
shallow = strcmp(config.equation,'shallow');

% % wave data and speed c
% uout  = uc(:, index);
% par.c = uc(end,index);
% xout  = x;

par.c = uout(end);

L = ceil(abs(xout(1)));         % domain size
N = length(xout);               % number of grid points
h = (2*L)/N;
uwave = uout(1:end-1);

if (strcmp(config.method, 'Chebyshev'))
    N = N + 2;
end

% % adjust c if we want (to standardize)

% par.c = 10;
% uout(end) = par.c;

symmetry_after = false;
% symmetry_after = true;

if symmetry_after
    config.symmetry = 'none';
else
    config.symmetry = 'L2squaredflip';
end

% adjust N or L if we want to
if exist('L_new','var')
    L = L_new;
end

% % if we change anything, need to send through fsolve again
% iter = 1000;
% [xout, uout] = fsolveequation(xout, uout, par, N, L, config, iter);
% uwave = uout(1:end-1);
% h = (2*L)/N;

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

% linear interpolation
Duhalf_fine = interp1(xhalf,Duhalf,xfine);
Duhalf_fine(isnan(Duhalf_fine)) = 0;

% zero derivative is where we have a sign change
zDer = find(diff(sign(Duhalf_fine)));
% but derivative is discontinuous, so we only only want every other one
numMinMax = 2;                      % how many min/max we want
zDer      = zDer(2*(1:numMinMax));  % take every other one, starting at second 
zDer_x    = xfine(zDer);            % x values of deriative

% % which min/max we want
% minmax = 1;
% % x value for join
% join_x = zDer_x(minmax);

% better way to do this, using known spacing
start = 1;
index = 2;

join_start = zDer_x(start);

join_x = join_start + (index - 1)*(spacing/2);

% add this line to find half-way waves
% join_x = (zDer_x(minmax) + zDer_x(minmax+1) )/ 2;

% % add this line with minmax = 1 to find pulse 0 (half way to first min)
% % we might not be able to do this for shallow water eq
% join_x = join_x / 2;

% where to join the waves
join_pt = round(join_x / h) - 1;

% right half-wave
ud_half    = uout( join_pt : center + join_pt - 1 );
ud_right   = flip( ud_half(1:end-1) );

% for periodic BCs, remove the final point
if strcmp(config.BC, 'periodic')
    ud_right   = ud_right(1:end-1);
end

ud         = [ ud_half ; ud_right; par.c ];

% if we use Chebyshev, interpolate back onto Chebyshev grid
if strcmp(config.method,'Chebyshev')
    ud = [interp1(x_uniform, ud(1:end-1), xout) ; par.c];
 
    % flip back around 
    xout = flip(xout);  
    ud = [ flip(ud(1:end-1)) ; par.c ];
end

% run joined pulse through Newton solver
% Newton solver on right half wave
iter = 100000;
[~, ud_out] = fsolveequation(xout, ud, par, N, L, config, iter);

% if we enforce symmetry afterwards, send back through Newton
% solver with symmetry enforced

if symmetry_after
    % save this in case we want to compare
    ud_old = ud_out
    % enforce symmetry
    config.symmetry = 'L2squaredflip';
    [~, ud_out] = fsolveequation(xout, ud_out, par, N, L, config, iter);
end

% plot double wave before and after Newton solver
figure;
plot(xout, ud(1:end-1), xout, ud_out(1:end-1));

legend('initial guess','Newton solver output');
title('double pulse');

findsymm(xout, ud_out, config);

x = xout;
u = ud_out;

end


