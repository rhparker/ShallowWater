% generates average eigenfunction from an eigenfunction
% with eigenvalue on the imaginary axis
%   average gives even function for real part
%   average gives odd function for imag part
% do this before fsolving

function vavg = avg_imag_eigenfn(v)

    vreal = real(vin);
    vimag = imag(vin);
    vreal_avg = [vreal(1); 0.5*(vreal(2:end) + flip(vreal(2:end)))];
    vimag_even = [vimag(1:N/2+1); -vimag(N/2+2:end)];
    vimag_even_avg = [vimag_even(1); 0.5*(vimag_even(2:end) + flip(vimag_even(2:end)))];
    vimag_odd_avg  = [vimag_even_avg(1:N/2+1); -vimag_even_avg(N/2+2:end)];
    vavg = vreal_avg + 1i * vimag_odd_avg;
end