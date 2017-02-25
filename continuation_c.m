% continuation code for wave speed c

% which equation to use
% use either shallow water equation or 5th order KdV
% config.equation = 'shallow';
config.equation = 'KdV';

% BCs to use
config.BC       = 'periodic';
% config.BC       = 'Neumann';

% which numerical method to use
% Fourier only works with periodic BCs
config.method   = 'Fourier';
% config.method   = 'fdiff';

% true-false parameters, for convenience
shallow = strcmp(config.equation,'shallow');
periodic = strcmp(config.BC,'periodic');
Fourier = strcmp(config.method,'Fourier');

% domain bounds and size of grid
% need more grid points for shallow water equation
if shallow
    L = 20;
    N = 2049;
else
%     L = 50;
    L = 100;
%     L = 200;
    N = 4097;             % Fourier
%     N = 1000;          % finite difference
end

% domain and step size
x = linspace(-L,L,N)';                    
h = 2*L / (N-1);

% deal with boundary conditions
if periodic
    % for periodic domain, remove final point
    % since we equate it with initial point
    x = x(1:end-1);
    N = N-1;
end

% generate differentiation matrices
if Fourier
    % use Fourier spectral methods for periodic domain
    % Fourier spectral method
    % generate differentiation matrices
    [D, D2, D3, D4, D5] = D_fourier(N, L);
else
    % finite difference method
    % generate differentation matrices
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);
end

% parameters and initial conditions
if shallow
    % parameter b to use; this will be fixed
    % for double pulses to exist, need b > 2 for this
    par.b = 2.5;

    % compute c based on b
    par.c = cstar(par.b);

    % we have exact solution for single pulse for these paramaters 
    % so start with that
    u = shallow_sol(x,par.b);
    
% 5th order KdV equation
else
    % we only have solution for one value of c, so we start there
    par.c = 36/169;
    u = (105/338)*sech( x / (2 * sqrt(13) ) ).^4;
end

% can use this as input for Newton solver
uin = [u; par.c];

%% some checks we can run

% % check to see that our pulse (numerically) solves the eq
[F,J] = integratedequation(u,par,N,config,D,D2,D3,D4,D5);

% % this should be 0 since u is a solution
% plot(x, F)

% % this should be 0 since D*u is an eigenvector with eigenvalue 0
% plot(x,J*D*u)

% % look at and plot eigenvalues around single pulse
% % eigenvalues distributed as expected, so that is good
% num = 50;
% center = -100;
% [lambda, V, J] = eigs_linear(x, u, b, par, L, N, num, center, 'integrated');

%% secant continuation code in parameter c

% number of iterations
iterations = 500;

% continuation parameters
contPar.numContSteps    = iterations;
contPar.Name            = 'c';  % continuation parameter continued in

if shallow
    contPar.ds          = 5;    % continuation step size: should always be positive!
    contPar.initial_ds  = 5;    % initial step: sign determines direction of continuation
else
    contPar.ds          = 2.0e-1;       % continuation step size: should always be positive!
    contPar.initial_ds  = 2.0e-1;       % initial step: sign determines direction of continuation
end

% system parameters
% initial condition is our exact solution from above
u0 = u;

%% find two initial points on the bifurcation curve

% first point
options = optimset('Display','iter','Algorithm','levenberg-marquardt','MaxIter',30,'Jacobian','on');
[u1,fval,exitflag,output,jacobian1]  = fsolve( @(u) integratedequation(u,par,N,config,D,D2,D3,D4,D5),u0,options);

v0 = [u1; getfield(par,contPar.Name)];  % v0 is the first point with parameter name

parameter = getfield(par,contPar.Name); % Set a few facility vectors
normu     = norm(u1);
contdata  = v0;

% second point
par = setfield(par,contPar.Name,getfield(par,contPar.Name)+contPar.initial_ds); % increase par by ds
[u2,fval,exitflag,output,jacobian1]  = fsolve(@(u) integratedequation(u,par,N,config,D,D2,D3,D4,D5),u1,options);

v1 = [u2; getfield(par,contPar.Name)]; % v1 is the second point with parameter name

parameter = [parameter getfield(par,contPar.Name)];
normu     = [normu     norm(u2)];
contdata  = [contdata  v1];


%% Continuation code
% At each continuation step
for i = 1:contPar.numContSteps

  % Predictor
  v = v1 + (v1 - v0)/norm(v1 - v0, 2) * contPar.ds;
  
  disp(['Predictor = ',num2str(v(end))]);
  % Call fsolve 
  [v,res,exitflag,output,jacobian2] = fsolve(@(v) FixedPointSecantPredictorCorrector(v,v1,v0,N,config,D,D2,D3,D4,D5,par,contPar),v,options); 
  
  disp(['Step = ',int2str(i),' Parameter = ',num2str(v(end)), ' Norm(residual,inf) = ',num2str(norm(res,inf))]);
  % Update output parameters
  parameter = [parameter v(end)];
%  normu     = [normu norm(v(1:end-1))];
  normu     = [normu norm(v(1:end-1))];
  contdata  = [contdata v];

  plot(parameter,normu,'-o');   % plot the bifurcation diagram on the fly
  xlabel(contPar.Name);
  ylabel('Norm of the solution');drawnow;
  
  % Prepare for the next continuation step
  v0 = v1;
  v1 = v;
end


