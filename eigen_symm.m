% construct symmetric eigenfunction by rotation in complex plane
% eigenvalue is (nearly) pure imaginary, so also
% use fsolve to remove this small imaginary part

function [ls, vs] = eigen_symm(l, v, x, J, config)
    % find best angle of rotation
    [vr,theta]   = rotate_evec(x, v, config);
    
%     % nothing really changes after just the rotation so don't need
%     max_rotate   = max( abs( J*vr - eVals(index)*vr));
%     max_real_flipdiff = max( real(vr(start:end)) - flip(real(vr(start:end))) );
%     max_imag_flipdiff = max( imag(vr(start:end)) + flip(imag(vr(start:end))) );

    [vs_real, vs_imag, ls] = eig_solve_symm(J, l, x, real(vr), imag(vr), config);
    ls = 1i * ls;
    vs = vs_real + 1i * vs_imag;
    
end