load 100F_double2
% load 100F;

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

%%

% cumulative integral of psi
% ditch first point since periodic and no equivalent
% on the R end
half_len = length(xout) / 2;
h = xout(2) - xout(1);
center = half_len + 1;
left_int = zeros(1, half_len);
right_int = zeros(1, half_len);

for index = 1:128
    left_int(index) = sum( psi_fd(1 : 1 + index) );
    right_int(index) = -sum( psi_fd( half_len + index : end) );
end

left_deriv = phi(1:center-1);
right_deriv = phi(center : end);

% Melnikov integral, before simplification
M1 = (left_int*left_deriv + right_int*right_deriv) * h;

% this is <q, q_c>
M2 = (uout(1:end-1)'*psi_fd) * h;

% this is q(0) * int q_c
M3 = uout(center) * sum(psi_fd) * h;
