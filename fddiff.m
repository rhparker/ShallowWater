%% generates finite difference matrix using Fornberg weights
%   N       : size of matrix
%   m       : order of derivative
%   h       : spatial step size
%   pts     : num of points to use on either side of center
%   BC      : boundary conditions ('periodic' for periodic)
% note that pts must be at least m/2 (rounded up) for this to work

function D = fddiff(N, m, h, pts, BC)
    wts = Fornberg_weights(0,-pts:pts,5);
    D   = sparse_multidiag(N, wts(m+1,:), BC)/(h^m);
end

% generates sparse multidiagonal matrix S
%   N      : size of matrix
%   vals   : value to put on the diagonals
%   config : if 'periodic' then make periodic matrix
% vals must have an odd number of elements; the middle 
% element goes on the main diagonal
function S = sparse_multidiag(N, vals, config)
    len = length(vals);
    % start with zero sparse matrix
    S = sparse(N,N);
    % if odd number of vals, keep going
    if mod(len,2) == 1
        diags = -(len-1)/2 : (len-1)/2;
        for i = 1:len
            S = S + sparse_diag(N, diags(i), vals(i), config);
        end
    end
end

% generates sparse diagonal matrix S
%   N      : size of matrix
%   d      : which diagonal to use (0 is main diagonal)
%   val    : value to put on the diagonal
%   config : if 'periodic' then make periodic matrix
function S = sparse_diag(N, d, val, config)
    diag = abs(d);
    S = sparse(1:N-diag,diag+1:N,val*ones(N-diag,1),N,N);
    % if we want a periodic matrix
    if exist('config','var')
        if (strcmp(config,'periodic')) && (diag ~= 0)
            S = S + sparse_diag(N, N-diag, val)';
        end
    end
    % if we specify a negative diagonal, take transpose
    if d < 0
        S = S';
    end
end
