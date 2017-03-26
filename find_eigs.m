%% eigenvalues of integrated operator

% single pulse
% load 0single_fourier;
% load 0single_fdiff;
% load 5single;

% load Sh4double1a; par.b=b;
load 2double2a;
uout = ud_out;
% load 6single;
% load 0singlefourier;
% load 0singleneumann; 

% default parameters
par.c = uout(end);              % wave speed
N = length(xout);               % current grid size
L = -xout(1);                   % current domain length;
h = 2*L/N;                      % current grid spacingn n b

% % if we want, we can change paramaters

% % change grid size
% N = 20000;                  % finite diff
% N = 2048;                   % Fourier
% N = 128;
% N = 4096;
% N = 400;

% % change domain size
% L = 100;

% % change speed c (to standardize between methods)
% par.c = 82.5;

% if we change stuff, run through Newton solver
if N ~= length(xout) || L ~= -xout(1) || par.c ~= uout(end)
    % interpolate onto a larger grid, or with a longer domain, or with
    % different c
    [xnew, unew] = fsolveequation(xout, uout, par, N, L, config, 10000);
else
    % if we don't interpolate
    xnew = xout; unew = uout;
end

% differentiation matrices
Fourier = strcmp(config.method,'Fourier');
if Fourier
    [D, D2, D3, D4, D5] = D_fourier(N, L);
else
    % grid spacing
    h = 2*L/N;
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC, 3);
end

% generate flipped and average wave
[uflip, uavg] = flip_avg_wave(unew, config);

% average pulse or not; this makes things nicer
average_pulse = true;
% average_pulse = false;

if average_pulse
    unew = uavg;
    % run the averaged pulse through fsolve
    [xnew, unew] = fsolveequation(xnew, unew, par, N, L, config, 10000);
end

% just the wave
uwave = unew(1:end-1);
[F,J] = equation(uwave,par,N,config,D,D2,D3,D4,D5);


% % afsolve with nonintegrated equation (5th order)
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

% [F,J] = equation(uavg(1:end-1),par,N,config,D,D2,D3,D4,D5);
% plot(xnew, F);

% [F,J] = integratedequation(uwave,par,N,config,D,D2,D3,D4,D5);
% plot(xnew, F)
% plot(xnew, J*D*uwave)

% % check eigenvalues of integrated equation
% num = 50;
% center = -100;
% [int_lambda, ~, ~] = eigs_linear(xout, uwave, par, config, num, center, 'integrated');
% plot(int_lambda, zeros(length(int_lambda)), '.');


%% eigenvalues

% eigenvalues of nonintegrated equation

% find the ideal exponential weight (at least I hope so)
% in terms of c
a = find_exp_wt(par.c);
% a = a / 2;
% a = 0.2;
a = 0;

% eigenvalue of the constant solution for linearization about zero solution
lambda_const = -a^5 + a^3 - par.c * a;

% % use eig
[lambda, V, J] = eig_linear(xnew, uwave, par, config, 'nonintegrated', a);

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
    cutoff = 0.1;
%     indices = find( abs(real(lambda)) > cutoff);
    indices = find( real(lambda) > cutoff);
    eVals = lambda(indices);
    eVecs = V(:, indices);
    integ = trapz(xnew,eVecs);
end

% imag_eval = false;
imag_eval = true;

if imag_eval
    % % grab the eigenvalue nearest the one we found from the 
    % % weighted space
    % target = 0.0215;
%     target = 0.6423;
%     target = 0.0149;     % 5double1a, N=256, c=4.4459
    target = 0.0015;     % 2double2a, N=256, c=9.4812
%     target = 0.0629;     % 2double1a, N=256, c=9.4812
%     target = 0.4596;   % 4double1a, N=256, c=32.6519
%     target = 0.2277;   % 3double1a, N=256, c=20.6361
    threshold = 0.001;
    index = 1;
    indices = find(abs(imag(lambda) - target) < threshold);
    eVals = lambda(indices);
    eVecs = V(:, indices);

%     % eliminate small real part with fsolve
%     [vout, lout] = eig_solve(J, i*imag(eVals(index)), eVecs(:,index), 'imag');
%     [vout, lout] = eig_solve(J, i*imag(eVals(index)), eVecs(:,index), 'fix_restrictnorm');

    % for now, we don't need to do that since we will
    % get rid of small real part when we fsolve for symmetry
    vout = eVecs; 
    lout = 1i*imag(eVals);   % take imaginary part only

    % average real(vout(x)) and real(vout(-x))
    vreal = real(vout);
    vreal_avg = [vreal(1); 0.5*(vreal(2:end) + flip(vreal(2:end)))];

    % start with the averaged real part
    vreal = vreal_avg;
    % compute the imaginary part from this
    vimag = (-1/imag(lout))*J*vreal;

    % fsolve and reconstruct eigenvector
    % fsolve here can change real part, imag part of eigenvector
    % fsolve here can change eigenvalue
    % at present, fsolve has no additional symmetry-enforcing conditions
    [v3_real, v3_imag, l3] = eig_solve_symm(J, lout, xnew, vreal, vimag, config);
    v3 = v3_real + 1i * v3_imag;
    max_before = max( abs( J*eVecs(:,index) - eVals(index)*eVecs(:,index)) );
%     max_after  = max( abs( J*vout - lout*vout));
    max_symm   = max( abs( J*v3 - 1i*l3*v3));

    max_real_flipdiff_before = max( real(eVecs(2:end)) - flip(real(eVecs(2:end))) );
    max_imag_flipdiff_before = max( imag(eVecs(2:end)) + flip(imag(eVecs(2:end))) );

    max_real_flipdiff = max( real(v3(2:end)) - flip(real(v3(2:end))) );
    max_imag_flipdiff = max( imag(v3(2:end)) + flip(imag(v3(2:end))) );
    
    imag_diff_before = max( imag(eVecs) - (-1/imag(eVals))*J*real(eVecs));
    imag_diff_symm   = max( imag(v3) - (-1/l3)*J*real(v3));
end
%% construct another eigenfunction
% 
% if imag_eval
%     [vout2, lout2] = eig_solve(J, i*imag(eVals(index)), vavg, 'fix_restrictnorm');
%     max_avg  = max( abs( J*vout2 - lout*vout2));    
% end

% % % check for symmetry
% if strcmp(config.BC, 'periodic')
%     eVecFlip = [ eVecs(1) ; flip(eVecs(2:end)) ];
% %     voutflip = [ vout(1) ; flip(vout(2:end)) ];
% end;
% 
% % plot(xnew, abs(eVecs) - abs(eVecFlip));
% max_flip = max( abs( abs(eVecs) - abs(eVecFlip) ) );

% % integrals
% integout = trapz(xnew,vout);
% integ = trapz(xnew,eVecs);

% figure;
% plot(xnew, exp(-a*xnew).*real(vout));
% title('real part of unweighted eigenfunction, eigenvalue 0.1534i');
% 
% figure;
% plot(xnew, exp(-a*xnew).*imag(vout));
% title('imag part of unweighted eigenfunction, eigenvalue 0.1534i');
% 
% figure;
% plot(xnew, imag(vout));
% title('imag part of eigenfunction, a = 0, eigenvalue 0.0215i');
% 
% figure;
% plot(xnew, real(vout));
% title('real part of eigenfunction, a = 0, eigenvalue 0.0215i');


% figure;
% plot(lambda, '.');
% plot_title   = ['Zoom of eigenvalues using eigs of double pulse 2, exp wt a = ',num2str(a)];
% method_title = [config.method,', ',config.BC,' (N = ',num2str(N),', L = ',num2str(L),') '];
% eig_title  = ['lambda = ',num2str(lambda(1))];
% title({plot_title, [method_title, ' ']}); 
% axis([-0.4 0.4 -3 3]);

% title('weighted space eigs for single pulse, known solution, fdiff/Neumann')
% axis([-0.002 0.003 -0.05 0.05]);

% figure;
% plot(lambda, '.');
% title('weighted space eigs for single pulse, known solution, fdiff/Neumann')
% axis([-0.001 0.001 -0.01 0.01]);

% figure;
% index = 3;
% plot(xnew,V(:,index))
% title(strcat('eigenfunction corresponding to eigenvalue:  ',num2str(lambda(index))))

% % use for complex eigenfunctions
% indices = 4:5;
% lambdaV = lambda(indices);
% legendcell = cellstr(num2str(lambdaV))
% 
% figure;
% plot( xnew, real( V(:,indices)) );
% legend(legendcell);
% title('Real part of eigenfunctions')
% 
% figure;
% plot( xnew, imag( V(:,indices)) );
% legendcell = cellstr(num2str(lambdaV))
% legend(legendcell);
% title('Imaginary part of eigenfunctions')
