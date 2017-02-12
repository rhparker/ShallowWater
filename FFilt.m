% Fourier filter
function F = FFilt(x, N, L)
    F = zeros(N);
    for row = 1:N
        for col = 1:N
            F(row, col) = g_j_sigma(x(row),x,col,N,L);
        end
    end
end

function g = g_j_sigma(x, x_grid, j, N, L)
    % constant c_n_sigma
    c_n_sigma = ones([1 N/2+1]);
    c_n_sigma(1) = 2; c_n_sigma(end) = 2;
    % n values we need
    n = 0:N/2;
    % scale x values to [0, 2 pi]
    z  = pi*(x/L + 1);
    zj = pi*(x_grid(j)/L + 1);
    % filtered values and cosine values for sum
    sigma   = filter( n / (N/2) );
    cosines = cos( (z - zj) * n );
    g = (2/N)*sum( 1./c_n_sigma .* sigma .* cosines);
end

% specific filter to use
function sigma = filter(nu)
    sigma = exp_filt(nu, 35, 2, 0.8);
end

% exponential filter
function sigma = exp_filt(nu, alpha, p, nu_crit)
    filtered_index = find( abs(nu) > nu_crit);
    filtered = zeros([1 length(nu)]);
    filtered(filtered_index) = ones([1 length(filtered_index)]);
    exp_filter = exp(-alpha* ( (abs(nu) - nu_crit)/(1 - nu_crit) ).^p);
    sigma = (1 - filtered) + filtered.*exp_filter;
end
