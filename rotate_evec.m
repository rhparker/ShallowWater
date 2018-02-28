% rotate eigenvector by unit complex number so that 
% real part is even; this should make imag part odd

function [u, theta] = rotate_evec(x, v, config)
    iter = 100;
    options=optimset('Display','iter','MaxIter',iter); 
    options.TolFun = 1e-10;
    options.TolX   = 1e-10;
    
    N = length(x);
    L = -x(1);
    h = 2*L/N;
    theta = 0;
    v_real = real(v);
    v_imag = imag(v);
    [output,fval] = fsolve( @(input) L2_flip_norm(x, v_real, v_imag, input, h, config), theta, options);
    theta = output;
    u_real = cos(theta/2).*v_real - sin(theta/2).*v_imag;
    u_imag = sin(theta/2).*v_real + cos(theta/2).*v_imag;
    u = u_real + 1i*u_imag;
end

% want real part to be even
% and imaginary part to be odd
function val = L2_flip_norm(x, v_real, v_imag, theta, h, config)
    N = length(x);
    if strcmp(config.BC, 'periodic')
        start = 2;
        center = N/2 + 1;
    else
        start = 1;
        center = (N+1)/2;
    end
    % real part
    u_r = cos(theta/2).*v_real - sin(theta/2).*v_imag; 
    % imaginary part
    u_i = sin(theta/2).*v_real + cos(theta/2).*v_imag;
    
    % deal with periodicity
    u_r = u_r(start:end);
    u_i = u_i(start:end);
    
    % check real part for even
    sum_even = sqrt( h * sum( (u_r - flip(u_r)).^2) );
    
    % check imag part for odd
    u_i_flip = [-u_i(1:center) ; u_i(center+1:end) ];
    sum_odd = sqrt( h * sum( (u_i_flip - flip(u_i_flip)).^2) );
    
    val = sum_even + sum_odd;
    val = sum_even;
end