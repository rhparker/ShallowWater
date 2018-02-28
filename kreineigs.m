% returns one of the Krein eigenvalues
% basically just calls the Krein matrix function
% useful for fsolve

function l = kreineigs(z, n, filename)
    [Kz, lambda] = kreinmatrix_z_noinv(z, filename);
    l = lambda(n);
end