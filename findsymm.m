% find symmetry of function
function symm_sum = findsymm(x, u, config)

% % use continuation data
% index = length(u(1,:));
% xout = x;
% uout = u(:,index);

% use regular data
uout = u;
xout = x;

c_start = 36/169;
uorig = (105/338)*sech( xout / (2 * sqrt(13) ) ).^4;
uorig = [uorig; c_start];

% paramaters
N = length(xout);
L = -xout(1);
h = (2*L)/N;

% assume we use Fourier, for now
[D, D2, D3, D4, D5] = D_fourier(N, L);

uwave = uout(1:end-1);
par.c = uout(end);
u2 = uout;

usum = uout(2:N);
flipdiff = abs(usum - flip(usum));
symm_sum  = h * sum( (usum - flip(usum)).^2 );
symm_sum2 = h * sum( (D*uwave).*(uwave - uorig(1:end-1)) );
F =  integratedequation(xout,uwave,par,N,config,D,D2,D3,D4,D5,uwave);

disp(['symmetry enforcment: ',config.symmetry]);
disp(['wave speed c: ',num2str(par.c)]);
disp(['integral of (u(x) - u(-x))^2: ',num2str(symm_sum)]);
disp(['integrated diff from symmetric fn: ',num2str(symm_sum2)]);
disp(['max of integrated operator: ',num2str(max(F(1:end-1)))]);
disp(['max of abs(u(x) - u(-x)): ',num2str(max(flipdiff))]);


% figure;
% uflip = flip_avg_wave(u2, config);
% plot(xout, u2(1:end-1) - uflip(1:end-1))
% plot(xout(2:end),usum - flip(usum));
% plot(xout(2:end), (usum - flip(usum)).^2 );

end