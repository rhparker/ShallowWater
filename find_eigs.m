%% eigenvalues of integrated operator

% single pulse
% load 0single;
% load 0single_fdiff;
% load 5single;

% load Sh4double1a; par.b=b;
% load 5double1;
% uout = ud_out;
% load 6single;
% load 0singlefourier;
% load 0singleneumann; 

load 100F;
index = 2;

uout = ud_out(:,index);

% uout = umc_2_4;

target = targets(index);
threshold = 0.0001;

% default parameters
par.c = uout(end);              % wave speed
N = length(xout);               % current grid size
L = ceil(abs(xout(1)));         % current domain length;
h = 2*L/N;                      % current grid spacing
N_old = length(xout);

% % if we want, we can change paramaters

% % change grid size
% N = 20000;                  % finite diff
% N = 2048;                   % Fourier
% N = 128;
% N = 4096;
% N = 400;

% % change domain size
% L = 100;

N = 400;

% % change speed c (to standardize between methods)
% par.c = 82.5;

config.form = 'integrated';

% if we use Chebyshev methods, add 2 to N since we have
% homogeneous Dirichlet BCs at endpoints
if strcmp(config.method,'Chebyshev')
    N = N+2;
    N_old = N_old + 2;
end

% config.symmetry = 'none';

% if we change stuff, run through Newton solver
if N ~= N_old || L ~= ceil(abs(xout(1))) || par.c ~= uout(end)
    % interpolate onto a larger grid, or with a longer domain, or with
    % different c
    if strcmp(config.method,'Chebyshev')
        [xnew, unew] = fsolveequation(xout, uout, par, N, L, config, 10000);
    else
        [xnew, unew] = fsolveequation(xout, uout, par, N, L, config, 10000);
    end
else
    % if we don't interpolate
    xnew = xout; unew = uout;
end

% differentiation matrices

if strcmp(config.method,'Fourier')
    [D, D2, D3, D4, D5] = D_fourier(N, L);
elseif strcmp(config.method,'Chebyshev')
    % for Chebshev methods we use Dirichlet BCs so 
    % we need to increase N by 2 to account for the
    % two boundary points
    [D, D2, D3, D4, D5] = D_cheb(N, L, config);
else
    % grid spacing
    h = 2*L/N;
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC, 3);
end

% average pulse or not; this makes things automatically symmetric
% average_pulse = true;
average_pulse = false;

if average_pulse
    % generate flipped and average wave
    [uflip, uavg] = flip_avg_wave(unew, config);
    unew = uavg;
    % run the averaged pulse through fsolve
    [xnew, unew] = fsolveequation(xnew, unew, par, N, L, config, 10000);
end

% just the wave
uwave = unew(1:end-1);

% need configuration without symmetry to find eigenvalues
config_nosymm = config;
config_nosymm.symmetry = 'none';

% % fsolve with nonintegrated equation (5th order)
% options = optimset('Display','iter','Algorithm','levenberg-marquardt','MaxIter',30,'Jacobian','on');
% [unew,fval,exitflag,output,jacobian1]  = fsolve( @(u) equation(u,par,N,config,D,D2,D3,D4,D5),uwave,options);
% unew  = [unew; par.c];

% % if we want to use the zero solution, uncomment this
% uwave = zeros([length(uwave) 1]);

%
% several checks that we can do
%

%
% test plots of difference between function and flipped function
% as well as derivative
%

% to test symmetry, we need to take this wave and flip it
% generate flipped wave u(-x)

% % plot of difference between wave and flipped wave
% figure;
% [uflip2, ~] = flip_avg_wave(unew, config);
% flipdiff = unew(1:end-1) - uflip2(1:end-1);
% plot(xnew, flipdiff);
% title('u(x) - u(-x), after averaging and fsolve, Fourier, N=256');

% % plot comparing difference to derivative (normalized)
% figure
% uderiv = D*uwave;
% plot(xnew, flipdiff/norm(flipdiff), xnew, uderiv/norm(uderiv))
% title('u(x) - u(-x) and derivative of single pulse (both normalized)');
% legend('u(x) - u(-x)', 'derivative of single pulse');

%
% % check to see that our pulse (numerically) solves the eq
% % and that derivative is eigenvector of J with eigenvalue 0
%

% [F,J] = equation(xnew, uwave,par,N,config_nosymm,D,D2,D3,D4,D5);
% plot(xnew, F);

% [F,J] = integratedequation(xnew,uwave,par,N,config,D,D2,D3,D4,D5);
% plot(xnew, F(1:end-1))

% plot(xnew, J*D*uwave)

% % check eigenvalues of integrated equation
% num = 50;
% center = -100;
% [int_lambda, ~, ~] = eigs_linear(xout, uwave, par, config, num, center, 'integrated');

% [int_lambda, V_int, J_int] = eig_linear(xnew, uwave, par, config_nosymm, 'integrated');
% pt_spec = int_lambda( find(int_lambda <= 1) );
% plot(pt_spec, zeros(length(pt_spec)), '.', 'MarkerSize', 10);
% plot(int_lambda, zeros(length(int_lambda)), '.');


%% eigenvalues

% eigenvalues of nonintegrated equation

% find the ideal exponential weight (at least I hope so)
% in terms of c
a = find_exp_wt(par.c);
% a = a / 2;
% a = 0.1;
a = 0;

% eigenvalue of the constant solution for linearization about zero solution
lambda_const = -a^5 + a^3 - par.c * a;

% % use eig
config_nosymm.form = 'nonintegrated';
[lambda, V, J] = eig_linear(xnew, uwave, par, config_nosymm, 'nonintegrated', a);

% use eigs
% num    = 5;
% center = 0.6423i;
% center = 0.0215i;
% [lambda, V, J] = eigs_linear(xnew, uwave, par, config, num, center, 'nonintegrated', a);

eig_plot = false;
% eig_plot = true;

if eig_plot
    figure;
    plot(lambda, '.');
    plot_title   = ['Eigenvalues using eigs of double pulse 1, exp wt a = ',num2str(a)];
    method_title = [config.method,', ',config.BC,' (N = ',num2str(N),', L = ',num2str(L),') '];
    % eig_title  = ['lambda = ',num2str(lambda(1))];
    axis([ -1 1 -1 1]);
    title({plot_title, [method_title, ' ']}); 
end

%% play with eigenvalues

% for exponentially weighted space, if we have a good
% separation, we can extract the eigenvalues which are
% to the R of the essential spectrum

if a ~ 0
    cutoff  = -1;
    indices = find(real(lambda) > cutoff);
    eVals = lambda(indices);
    eVecs = V(:, indices);
    
    % for eigenvalues close to imag axis, we can fsolve them to there
    index = 2;
%     [vout, lout] = eig_solve(J, i*imag(eVals(index)), eVecs(:,index), 'fix_restrictnorm');
%     max_before = max( abs( J*eVecs(:,index) - eVals(index)*eVecs(:,index)) );
%     max_after  = max( abs( J*vout - lout*vout));
    

% otherwise grab the eigenvalues off the real axis
else
    % looking for real eigenvalues
    if imag(target) == 0
        cutoff = 0.0001;
        indices = find( abs( abs(real(lambda)) - target ) < cutoff);
        eVals = lambda(indices);
        eVecs = V(:, indices);
        integ = trapz(xnew,eVecs);
        
    % otherwise looking for imaginary eigenvalues
    else
        % % grab the eigenvalue nearest the one we found from the 
        % % weighted space
    
        % this finds both complex conjugates
        indices = find (abs( abs(imag(lambda)) - imag(target) ) < threshold);

%         this just finds one with pos imag part
%         indices = find (abs( imag(lambda) - imag(target) ) < threshold);

        eVals = lambda(indices);
        eVecs = V(:, indices);

    %     % eliminate small real part with fsolve
    %     [vout, lout] = eig_solve(J, i*imag(eVals(index)), eVecs(:,index), 'imag');
    %     [vout, lout] = eig_solve(J, i*imag(eVals(index)), eVecs(:,index), 'fix_restrictnorm');
    %     max_after  = max( abs( J*vout - lout*vout));
    
        % want eigenvalue with positive imag part, for consistency
        index = 1;

        % for now, we don't need to do that since we will
        % get rid of small real part when we fsolve for symmetry
        vout = eVecs(:,index); 
        lout = 1i*imag(eVals);   % take imaginary part only
        
        if strcmp(config.BC, 'periodic')
            start = 2;
        else
            start = 1;
        end
        
        % symmetry and max of eigenvalue problem before doing anything
        max_before = max( abs( J*eVecs(:,index) - eVals(index)*eVecs(:,index)) );
        max_real_flipdiff_before = max( real(eVecs(start:end)) - flip(real(eVecs(start:end))) );
        max_imag_flipdiff_before = max( imag(eVecs(start:end)) + flip(imag(eVecs(start:end))) );
        imag_diff_before = max( imag(eVecs(:,index)) - (-1/imag(eVals(index)))*J*real(eVecs(:,index)));

    %     % start with averaged real part and construct symmetric eigenvector
    %     % average real(vout(x)) and real(vout(-x))
    %     vreal = real(vout);
    %     vreal_avg = [vreal(1); 0.5*(vreal(2:end) + flip(vreal(2:end)))];
    %     % start with the averaged real part
    %     vreal = vreal_avg;
    %     % compute the imaginary part from this
    %     vimag = (-1/imag(lout))*J*vreal;
    % 
    %     % fsolve and reconstruct eigenvector
    %     % fsolve here can change real part, imag part of eigenvector
    %     % fsolve here can change eigenvalue
    %     % at present, fsolve has no additional symmetry-enforcing conditions
    %     [v3_real, v3_imag, l3] = eig_solve_symm(J, lout, xnew, vreal, vimag, config);
    %     v3 = v3_real + 1i * v3_imag;
    %     max_symm   = max( abs( J*v3 - 1i*l3*v3));
    %     max_real_flipdiff = max( real(v3(2:end)) - flip(real(v3(2:end))) );
    %     max_imag_flipdiff = max( imag(v3(2:end)) + flip(imag(v3(2:end))) );
    %     imag_diff_symm   = max( imag(v3) - (-1/l3)*J*real(v3));
    
        % rotate and get rid of imaginary part
    	[l5, v5] = eigen_symm(eVals(index), eVecs(:,index), xout, J, config);
   
        % symmetry and max info after
        max_rotate_fsolve   = max( abs( J*v5 - l5*v5));
        max_real_flipdiff_fsolve = max( real(v5(start:end)) - flip(real(v5(start:end))) );
        max_imag_flipdiff_fsolve = max( imag(v5(start:end)) + flip(imag(v5(start:end))) );
    end
end

