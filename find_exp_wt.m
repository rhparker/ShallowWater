% find optimal exponential weight for 5th order KdV
%   specify paramater c which is wave speeed
%   method specifies how we measure distance between eigenvalues
function alpha = find_exp_wt(c, method)
    if ~exist('method','var')
        method = 'real'
    end
    % use fsolve, starting at lambda = 0
    best_lambda = fsolve( @(lambda) top_roots_diff(c,lambda,method), 0);
    % find roots corresponding to this lambda
    nu = roots([ 1 0 -1 0 c -best_lambda ]);
    % third root is the one we want
    % should be sorted, but consider sorting it
    alpha = real(nu(3));
end

% difference between 3rd and 4th highest real parts
% of roots of char poly
% method = 'real' : distance between real parts
% method = 'abs'  : absolute value of distance
function d = top_roots_diff(c, lambda, method);
    r = sort_roots(c, lambda);
    % for absolute value, take dist between 3rd and 4th 
    % root, sorted by real part
    if strcmp('method', 'abs')
        d = abs( r(4) - r(3) );
    % otherwise, only need to look at real parts
    else
        % take unique real parts, and subtract last two
        r_real = unique(real(r));
        d = r_real(end) - r_real(end-1);
    end
end

% returns roots sorted by real part
function r = sort_roots(c, lambda);
    nu     = poly_roots(c, lambda);
    nu_sep = sortrows([real(nu) imag(nu)]);
    r      = to_complex(nu_sep);
end

% convert 2-tuple vectors to complex numbers
function c = to_complex(v)
    c = v(:,1) + v(:,2)*i;
end

% find roots of characteristic polynomial
function nu = poly_roots(c, lambda)
    nu = roots([ 1 0 -1 0 c -lambda ]);
end