% Krein matrix for Kdv5
% new version using Bloch wave formulation

%% setup
clear all;

% load wave data
load 100F;

config.symmetry = 'none';

% which solution to linearize about
% u = uout;
u = ud_out_2;
% u = um_2_4;

% u = uout;
% neg_directions = 2;
% pulses = 1;

% Bloch wave parameter / exponential weight
% mu = 0.006125 give us Krein bubble
% mu = 0.01           % good for Bloch wave
mu = 0.01;         % good for exp wt
% mu = 0.006125;

% Bloch wave parameter (this is imaginary)
% config.weight = -1i * mu;

% exponential weight (this is real)
config.weight = -mu;

% values of z to use for Krein matrix calculation
% range of z values
% zm = -0.21; zp = 0.21;
zm = 0;
zp = 0.12;
Nz = 100;
z = linspace(zm,zp,Nz + 1);         % interval over which to evaluation matrix
deltaz = z(2) - z(1);
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

% % take subspace S to be initial vectors in VA0
threshold = 0.1;
indices = find(abs(evalsA0) < threshold);
S = [ VA0(:, indices) ] ;
nu = evalsA0(indices);
num = length(indices);

% projection matrices
% easiest way to get a basis is to tack on standard
% basis vectors and use Gram-Schmidt via QR

Id = eye(N);

% M is change of basis matrix, which is orthogonal
M = [S Id(:, num:end-1)];
[M, ~] = qr(M, 0);

% PM is projection matrix in M-coordinates
% gets rid of v1, v2, Ker D
PMperp = Id; 
PMperp(:,1:4) = zeros(N, 4);
PM = Id - PMperp;

% PS is projection matrix in in std coord
PS = M*PM*M';
PSperp = M*PMperp*M';

% data for Krein eigenvalues
Kreineval = zeros(num, zsize);
Kreindets = zeros(1, zsize);
Kreindiagdets = zeros(1, zsize);
Kreinmatrices = zeros(num, num, zsize);
Krein1 = zeros(num, num, zsize);
Krein2 = zeros(num, num, zsize);
Kreinapprox = zeros(num, num, zsize);
Kreinapproxeval = zeros(num, zsize);
% Kreinderivs = zeros(num, num, zsize);
KS2norms = zeros(1, zsize);
KS1norms = zeros(1, zsize);
KSnorms = zeros(1, zsize);

% approximation to PSperp P(lambda) PSperp using A0
PA0M = M'*A0*M;
PA0M = blkdiag(Id(1:num, 1:num), PA0M(num+1:end, num+1:end));
PA0Minv = inv(PA0M);
PA0Minv = blkdiag(zeros(num), PA0Minv(num+1:end, num+1:end));
PA0inv = M*PA0Minv*M';

% Krein matrix
for kk=1:zsize
    Pz = A0 + 1i*z(kk)*A1;
    Ks1 = S'*(Pz*S);
    
    P1 = PSperp*Pz;
    P2 = PSperp*Pz*PSperp;
    
    % find inverse of P2; effectively we are doing this on Sperp
    P2M = M'*Pz*M;
    P2M = blkdiag(Id(1:num, 1:num), P2M(num+1:end, num+1:end));
    
    % invert with inv
    P2Minv = inv(P2M);
    P2Minv = blkdiag(zeros(num), P2Minv(num+1:end, num+1:end));
    P2inv = M*P2Minv*M';
    
%     % can use this version if you assume lambda = i z
%     Ks2 = (P1*S)'*(P2inv*P1*S);

    % this version works for any lambda
    Ks2 = S' * (Pz*PSperp*P2inv*PSperp*Pz*S);
    
    % approximate version
    Ks2approx = (P1*S)'*(PA0inv*P1*S);
    
% %     solve with linsolve
%     yM = linsolve(P2M, PMperp*M*Pz*S);
%     Ks2 = (P1*S)'*(M'*yM);

    % Krein matrix construction
    Ks = Ks1 - Ks2;
%     Ks = Ks1;

%     % derivative of Krein matrix
%     RR = Pz * PSperp * P2inv * PSperp * (1i * A1);
%     SS = Pz * PSperp * P2inv * (1i * A1) * P2inv * PSperp * Pz;
%     DKs = S' * (1i*A1 + RR + RR' - SS)*S;
%     Kreinderivs(:,:,kk) = DKs;

    % find Krein eigenvalues and determinant
    Kreineval(:,kk) = eig(Ks);
    Kreinmatrices(:,:,kk) = Ks;
    Krein1(:,:,kk) = Ks1;
    Krein2(:,:,kk) = Ks2;
    Kreinapprox(:,:,kk) = Ks1 - Ks2approx;
    Kreinapproxeval(:,kk) = eig(Ks1 - Ks2approx);
    Kreindets(1,kk) = det(Ks);
    Kreindiagdets(1,kk) = det(diag(diag(Ks)));
    KS2norms(1,kk) = norm(Ks2);
    KS1norms(1,kk) = norm(Ks1);
    KSnorms(1,kk) = norm(Ks);
end

% plot Krein eigenvalues
for kk=1:num
    figure;
    plot(z, real(Kreineval(kk,:)), '.', z, real(Kreinapproxeval(kk,:)), '.');
end

%% polynomials for determinant of Krein matrix
% % involve varying degrees of approximation
% 

% % enforce A1 to be skew symmetric
% % only consider doing for block wave formulation
% A1 = (A1 - A1')/2;

% first polynomial; set up for Bloch wave only
alpha11 = S(:,1)'*A1*S(:,1);
alpha12 = S(:,1)'*A1*S(:,2);
alpha22 = S(:,2)'*A1*S(:,2);
talpha11 = 1i*alpha11;
talpha22 = 1i*alpha22;

nu2 = nu(2);

beta11 = (PSperp*A1*S(:,1))'*(PA0inv*PSperp*A1*S(:,1));
beta12 = (PSperp*A1*S(:,1))'*(PA0inv*PSperp*A1*S(:,2));
beta22 = (PSperp*A1*S(:,2))'*(PA0inv*PSperp*A1*S(:,2));

beta11 = (PSperp*A1*S(:,1))'*(PA0inv*PSperp*A1*S(:,1));
beta12 = (PSperp*A1*S(:,1))'*(PA0inv*PSperp*A1*S(:,2));
beta22 = (PSperp*A1*S(:,2))'*(PA0inv*PSperp*A1*S(:,2));

temp = S'*(A0*PSperp*PA0inv*PSperp*A1*S);

c1 = talpha11*nu2;
c2 = -alpha12*conj(alpha12) + talpha11*talpha22 - beta11*nu2;
c3 = -talpha22*beta11 - 1i*conj(alpha12)*beta12 + 1i*alpha12*conj(beta12) - talpha11*beta22;
c4 = -beta12*conj(beta12) + beta11*beta22;

p = roots([c4 c3 c2 c1 0]);

% simpler polynomial

c0 = -abs(alpha12)^2 - beta11*nu2;
c1 = 2 * imag( conj(alpha12)*beta12 );
c2 = beta11*beta22 - abs(beta12)^2;

p2 = roots([c2 c1 c0]);

c0 = -nu2;
c1 = 0;
c2 = beta22;

p3 = roots([c2 c1 c0]);

r = zeros(1, num);
for kk = 1:num
    b = (PSperp*A1*S(:,kk))'*(PA0inv*PSperp*A1*S(:,kk));
    r(kk) = sqrt( nu(kk) / b );
end


