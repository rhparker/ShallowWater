%% shallow water wave equation
% Sanstede (2015), Chapneys and Groves (1997)

% BCs to use
config.BC       = 'periodic';

% domain bounds and size of grid
L = 20;
N = 2049;

% domain and step size
x=linspace(-L,L,N)';                    
h=2*L / (N-1);

% for periodic domain, remove final point
% since we equate it with initial point
x = x(1:end-1);
N = N-1;

% Fourier spectral method
% generate differentiation matrices
[D, D2, D3, D4, D5] = D_fourier(N, L);

% % finite difference method
% % generate differentation matrices
% [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);

% parameter b to use; this will be fixed
% for double pulses to exist, need b > 2 for this
b = 2.5;
b = 2.05;

% % this gives you approximately c = 1
% b = 2.29875;

% compute c based on b
par.c = cstar(b);

% we have exact solution for single pulse for these paramaters 
% so start with that
u = shallow_sol(x,b);

% can use this as input for Newton solver
uin = [u; par.c];

%% some checks we can run

% % check to see that our pulse (numerically) solves the eq
[F,J] = integratedshallow(u,b,par,N,D,D2,D3,D4);

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
% b will remain fixed

% continuation parameters
contPar.numContSteps    = 10;
contPar.ds              = 1;    % continuation step size: should always be positive!
contPar.initial_ds      = 1;    % initial step: sign determines direction of continuation
contPar.Name            = 'c';       % continuation parameter continued in

contPar.numContSteps    = 50;
contPar.ds              = 0.2;    % continuation step size: should always be positive!
contPar.initial_ds      = 0.2;    % 

% system parameters
% initial condition is our exact solution from above
u0 = u;

%% find two initial points on the bifurcation curve

% first point
options = optimset('Display','iter','Algorithm','levenberg-marquardt','MaxIter',30,'Jacobian','on');
[u1,fval,exitflag,output,jacobian1]  = fsolve( @(u) integratedshallow(u,b,par,N,D,D2,D3,D4),u0,options);

v0 = [u1; getfield(par,contPar.Name)];  % v0 is the first point with parameter name

parameter = getfield(par,contPar.Name); % Set a few facility vectors
normu     = norm(u1);
contdata  = v0;

% second point
par = setfield(par,contPar.Name,getfield(par,contPar.Name)+contPar.initial_ds); % increase par by ds
[u2,fval,exitflag,output,jacobian1]  = fsolve(@(u) integratedshallow(u,b,par,N,D,D2,D3,D4),u1,options);

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
  [v,res,exitflag,output,jacobian2] = fsolve(@(v) FixedPointSecantPredictorCorrector(v,v1,v0,b,L,N,D,D2,D3,D4,par,contPar),v,options); 
  
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



