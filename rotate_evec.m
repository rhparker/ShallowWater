% rotate eigenvector by unit complex number so that 
% real part is even; this should make imag part odd

function [u, theta] = rotate_evec(x, v, config)
    iter = 100;
    options=optimset('Display','iter','MaxIter',iter); 
    options.TolFun = 1e-20;
    options.TolX   = 1e-20;
    
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

function val = L2_flip_norm(x, v_real, v_imag, theta, h, config)
    N = length(x);
    if strcmp(config.BC, 'periodic')
        start = 2;
    else
        start = 1;
    end
    v = cos(theta/2).*v_real - sin(theta/2).*v_imag; 
    vsum = v(start:end);
    val = sqrt( h * sum( (vsum - flip(vsum)).^2) );
end