% use fsolve to find eigenvectors and eigenvalues
% use this in special case where eigenvalue is on imag axis
%  real part of eigenfunction is even
%  imag part of eigenfunction is odd
%  config is our standard config file (need to know BCs)
%  v_real is real part of eigenfunction

function [vreal, vimag, lout] = eig_solve_symm(J, l, x, v_real, v_imag, config)
    iter = 100;
    options=optimset('Display','iter','MaxIter',iter); 
    options.TolFun = 1e-20;
    options.TolX   = 1e-20;
    
    input = [v_real; v_imag; imag(l) ];
    [output,fval] = fsolve( @(input) eig_eq_symm(J, l, x, input, config), input, options);
    
    N = length(x);
    vreal = output(1:N);
    vimag = output(N+1:2*N);
    lout = output(end);
end

% eigenvalue equation
%    even extension for real part
%    odd extension for imaginary part
function f = eig_eq_symm(J, l, x, input, config)
%     % this part does even extension, remove for now
%     N = length(J);
%     if strcmp(config.BC, 'periodic')
%         start = 2;
%     else
%         start = 1;
%     end
%     % even extension for real part
%     uhalf_real = real(uhalf_real);
%     u_real = [uhalf_real; flip(uhalf_real(start:end-1))];

    N = length(x);
    L = -x(1);
    h = 2*L/N;

    % extract stuff from input
    u_real = input(1:N);
    u_imag = input(N+1:2*N);
    beta = input(end);

%     % we know how to compute the imaginary part
%     u_imag = -(1/beta)*J*u_real;

    u = u_real + i*u_imag;
    % to enforce symmetry take discrete L2 norm of
    % u_real(x) - u_real(-x)
    u_real_sum = u_real(2:end);
    symm_sum = sum( (u_real_sum - flip(u_real_sum)).^2 );
%     f = [J*u - l*u];

    % separate real and imaginary parts, so everything is real
%     f = [J*u_real + beta*u_imag; J*u_imag - beta*u_real; h*symm_sum];
    f = [J*u_real + beta*u_imag; J*u_imag - beta*u_real];
end
