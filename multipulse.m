%% generate multiple pulses
% this uses the already constructed double pulse

% load from continuation data

% load pregenerated single/double pulse

load 100F;

% which double pulse do we want
uout = ud_out_1;

% how many pulses we want
num_waves = 4;

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

% adjust N or L if we want to

% N = 1024;
% L = 50;
% h=2*L/N;

symmetry_after = false;
% symmetry_after = true;

if symmetry_after
    config.symmetry = 'none';
else
    config.symmetry = 'L2squaredflip';
end

% % if we change anything, need to send through fsolve again
% [xout, uout] = fsolveequation(x, uout, par, N, L, config);
% uwave = uout(1:end-1);
% h = (2*L)/N;

% chebychev points listed backwards, so flip around 
% for double pulse construction
if strcmp(config.method,'Chebyshev')
    xout = flip(xout);
    uout = [ flip(uout(1:end-1)) ; par.c ];
end

%% make multi-pulse

% where to split depends on whether we have periodic BCs or not
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
    uout = [spline(xout, uout(1:end-1), x_uniform) ; par.c];
else
    x_uniform = xout;
end

% interpolate onto a finer grid so we actually have peaks
N_fine = 50001;
center_fine = (N_fine + 1) / 2;
x_fine = linspace(-L, L, N_fine);
u_fine = spline(x_uniform, uout(1:end-1), x_fine);

% find highest peaks (peaks of double pulses)
[pks, locs] = findpeaks(u_fine);
[~, index] = sort(pks);
locs = sort( locs(index(end-1:end)) );

% left and right half waves
u_left    = u_fine( 1 : center_fine )';
u_right   = u_fine( center_fine + 1 : end )';

u_middle  = [u_fine( center_fine+1 : locs(2) )  u_fine( locs(1)+1 : center_fine) ]';

um_center = repmat(u_middle, num_waves-2, 1);

% construct multiple pulses
um = [ u_left ; um_center; u_right ];

% output too long, so need to chop off the ends
N_m = length(um);
center_m = (N_m + 1)/2;
center_diff = center_m - center_fine;
um = um( center_diff + 1 : center_diff + N_fine );

% interpolate back onto original grid
um = interp1(x_fine, um, xout);

% if we used Chebyshev, flip everything back around
if strcmp(config.method,'Chebyshev')
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
title('multipulse solution');

findsymm(xout, um_out, config);


