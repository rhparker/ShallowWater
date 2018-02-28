function Q = KDVsteady(q)  %5th-order KdV steady-state problem

global L N x b c ks
global D1 D2 D4

Qt = 2/15*ks^2*D4*q - b*ks*D2*q - c*q + 3*q.^2/2 ...
    + ks/2*(D1*q).^2 + ks*q.*(D2*q);
Q = real(Qt(1:N));

