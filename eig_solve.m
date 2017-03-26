% use fsolve to find eigenvectors and eigenvalues

function [vout, lout] = eig_solve(J, l, vin, config)
    if ~exist('config','var')
        config = 'fix';
    end

    iter = 100;
    options=optimset('Display','iter','MaxIter',iter); 
    options.TolFun = 1e-20;
    options.TolX   = 1e-20;
    
    % fix the eigenvalue, i.e. fsolve cannot change it
    if strcmp(config, 'fix')
        [vout,fval] = fsolve( @(u) eig_eq(J, l, u), vin, options);
        lout = l;

    % fix the eigenvalue and restrict norm to 1
    elseif strcmp(config, 'fix_restrictnorm')
        [vout,fval] = fsolve( @(u) eig_eq_restrictnorm(J, l, u), vin, options);
        lout = l;

    elseif strcmp(config, 'symm')
        [vout,fval] = fsolve( @(u) eig_eq_symm(J, l, u), vin, options);
        lout = l;

    % restrict eigenvalue to imaginary axis, but
    % allow fsolve to change it; for this we assume
    % the starting eigenvalue l is on or almost on
    % the imaginary axis and split eigenvalue problem
    % into real and imag parts
    elseif strcmp(config, 'imag')
        input = [real(vin); imag(vin); imag(l) ]
        [output,fval] = fsolve( @(u) eig_imag_axis(J, u), input, options);
        N = length(J);
        vout = output(1:N) + i*output(N+1:2*N);
        lout = i*output(end);
    end
end

% eigenvalue equation
function f = eig_eq(J, l, u)
    f = J*u - l*u;
end

% eigenvalue equation, restrict norm to 1
function f = eig_eq_restrictnorm(J, l, u)
    f = [J*u - l*u; norm(u) - 1.0];
end

% variant of eigenvalue equation, restrict eigenvalue to imaginary axis
% this is separate eigenvalue equations for real and
% imaginary part
function f = eig_imag_axis(J, input)
    N = length(J);
    u_real = input(1:N);
    u_imag = input(N+1:2*N);
    u = u_real + i*u_imag;
    l = input(end);
%     f = ((J*u_real + l*u_imag).^2 + (J*u_imag - l*u_real).^2).^(1/2);
    f = [J*u_real + l*u_imag; J*u_imag - l*u_real]
end

