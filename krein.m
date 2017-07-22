% linearization about 

load 100F_double2

% find the kernel of the 4th order operator
% 4th order operator is called L in Kapitula
% L is matrix J_int here
threshold = 1e-10;
index = find(abs(lambda_int) < 1e-10);
phi = V_int(:,index);

% we should be able to solve equation
% JL psi = phi
% where JL is 5th order operator as in Kapitula
% JL is matrix J here

% we also know that phi is the derivative of the stationary
% solution, so put that in for phi instead of the eigenfunction
% we found with eig
phi = D*uout(1:end-1);

% solve for psi using linsolve
psi = linsolve(D*J_int, phi);

% but we also know that psi is -du/dc, so we should be 
% able to compute that with finite differences
par.c = uout(end);
N = length(xout);
L = abs(xout(1));
iter = 1000;

step_size = 0.001;
par.c = uout(end) + step_size;
[~, uout_r] = fsolveequation(xout, uout, par, N, L, config, iter);

par.c = uout(end) - step_size;
[~, uout_l] = fsolveequation(xout, uout, par, N, L, config, iter);

psi_fd = -( uout_r(1:end-1) - uout_l(1:end-1) ) / (2 * step_size);

plot(xout, psi_fd);

% compare how well these two work

% matrix D in Kapitula is 1x1 here, so all we 
% have to do is compute <psi, L psi>
% L psi = u*, our stationary solution

d = psi'*uout(1:end-1);

d_fd = psi_fd'*uout(1:end-1);