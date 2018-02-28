% Krein matrix for Kdv5
% new version using Bloch wave formulation

%% setup
clear all;

% load wave data
load 100F;

config.symmetry = 'none';

% which solution to linearize about
u = ud_out_2;
neg_directions = 4;
pulses = 2;

% u = uout;
% neg_directions = 2;
% pulses = 1;

% Bloch wave parameter
% mu = 0.006125 give us Krein bubble
mu = 0.01;
% mu = 0.006125;
config.weight = -1i * mu;

% values of z to use for Krein matrix calculation
% range of z values
% zm = -0.21; zp = 0.21;
zm = 0;
zp = 0.12;
Nz = 500;
z = linspace(zm,zp,Nz + 1);         % interval over which to evaluation matrix
zsize = length(z);

% extract parameters from data
par.c = u(end);                 % wave speed
u = u(1:end-1);                 % just the wave
N = length(xout);               % current grid size
L = ceil(abs(xout(1)));         % current domain length;
h = 2*L/N;                      % current grid spacing


% generate differentiation matrix for Block wave formulation
[D, D2, D3, D4, D5] = D_fourier(N, L);
Dmu = D - config.weight*eye(N);

% matrices A1 and A0 needed for Krein matrix construction
% A1 is inverse of diff op
A1 = -inv(Dmu);

% A0 is (symmetric) linear operator representing energy
[F,A0] = integratedKdV(xout,u,par,config,D,D2,D3,D4,D5,u);

% eigenvalues of 5th order operator, for comparison
[~, JL] = KdV(xout,u,par,config,D,D2,D3,D4,D5,config.weight);
[JLev, JLeigs] = eig(JL);
JLeigs = diag(JLeigs);

% eigenvalues of A0
% should come sorted
[VA0,evalsA0] = eig(A0);
evalsA0 = diag(evalsA0);

%% Krein matrix construction

% if we want to project onto only the small eigenvalues
% we need to move some things around
% this swaps the first and second group of eigenstuff
% so that small ones are first, negative ones second
ind = pulses + 1;
num = pulses - 1;
evalsA0 = [ evalsA0(ind:ind+num) ; evalsA0(1:pulses) ; evalsA0(ind+pulses:end) ];
VA0 = [ VA0(:, ind:ind+num)  VA0(:, 1:pulses)  VA0(:, ind+pulses:end) ];


% block matrices
% numA0 = neg_directions;
numA0 = 1;

evalsA0abs = abs(evalsA0);
D0m12 = diag(sqrt(1./evalsA0abs));
A1new = D0m12*VA0'*A1*VA0*D0m12;
A0m = -eye(numA0);
A0p = eye(N-numA0);
A1m = A1new(1:numA0,1:numA0);
A1p = A1new(numA0+1:N,numA0+1:N);
A1L = A1new(numA0+1:N,1:numA0);

% Krein eigenvalues

Kreineval = zeros(numA0, zsize);
Kreindets = zeros(1, zsize);
Kreinmatrices = zeros(numA0, numA0, zsize);

for kk=1:zsize
    Ks = zeros(numA0);  
    Ks = (A0p + z(kk)*1i*A1p)\A1L;
    Ks = A1L'*Ks;
    Ks = A0m + z(kk)*1i*A1m  - z(kk)^2*Ks;
    % don't do this for now
%     Ks = -z(kk)*(Ks' + Ks)/2;
    Kreineval(:,kk) = eig(Ks);
    Kreinmatrices(:,:,kk) = Ks;
    Kreindets(1,kk) = det(Ks);
end  

