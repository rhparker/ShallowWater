function [dist, ht_diffs] = peak_dist_finterp(x, data)

    dist = [];
    ht_diffs = [];
%     vels = [];
    
    data_size = size(data);
    data_length = data_size(2);
    
    N = length(x);                  % current grid size
    L = -x(1);                      % current domain length
    h = 2*L/N;                      % current grid spacing
    
    % fourier frequencies
    delta_k = 2*pi/(N*h);
    k = [0:delta_k:N/2*delta_k,-(N/2-1)*delta_k:delta_k:-delta_k]';
    
    
    for index = 1:data_length
        % first find out roughly where the two peaks are        
        uout = data(:,index);
        [pks, locs] = findpeaks(uout);
        [d1, d2] = sort(pks);
        pklocs = locs(d2(end-1:end));
        % left and right peak locations
        xL = x(min(pklocs));
        xR = x(max(pklocs));
        
        % FFT of derivative of function
        uhat = ifft(uout);
        uhatderiv  = 1i*k.*uhat;
        uhat2deriv = -k.^2.*uhat;
        
        options=optimset('Display','none');
%         options.TolFun = 1e-16;
%         options.TolX = 1e-16;
        % left zero
        zL = fsolve( @(x) tpoly(x, k, L, uhatderiv), xL, options);
        % right zero
        zR = fsolve( @(x) tpoly(x, k, L, uhatderiv), xR, options);
        % peak distance
        dist = [dist (zR-zL)];
        
        % difference in peak heights
        ht_diff = tpoly(zR, k, L, uhat) - tpoly(zL, k, L, uhat);
        ht_diffs = [ht_diffs ht_diff];
        
        
    end
end

function f = tpoly(x, k, L, coeffs)
    f = real( sum( coeffs.*exp(-1i*k*(x+L)) ) );
end
