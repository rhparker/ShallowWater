% use fsolve to find zeros of krein matrix
% use the noinvert version since faster

filename = '1500F_double2_1195_krein';
% filename = '100F_double2_krein';

l = 0.0050 + 3.9867i;
% l = 0.0691i
% l = 1.5;
start = -l^2;
start = [real(start) imag(start)];

options = optimset('Display','iter');
options.TolFun = 1e-15;
options.TolX = 1e-15;

[z, fval] = fsolve( @(z) abs(det(kreinmatrix_z_noinv(z(1) + z(2)*1i, filename))), start, options);
% [z, fval] = fsolve( @(z) kreineigs(z, 1, filename), start, options);

[k, l] = kreinmatrix_z_noinv(z(1) + z(2)*1i, filename);
[kr, lr] = kreinmatrix_z_noinv(z(1), filename);