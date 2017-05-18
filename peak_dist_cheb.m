function [dist, ht_diffs] = peak_dist_cheb(x, data)

    dist = [];
    ht_diffs = [];
    
    data_size = size(data);
    data_length = data_size(2);
    
    N = length(x);                  % current grid size
    L = ceil(abs(x(1)));            % current domain length
    h = 2*L/N;                      % current grid spacing
    
    for index = 1:data_length
        % first find out roughly where the two peaks are        
        uout = data(:,index);
        [pks, locs] = findpeaks(uout);
        [d1, d2] = sort(pks);
        pklocs = locs(d2(end-1:end));
        % left and right peak locations
        xL = x(min(pklocs));
        xR = x(max(pklocs));
        
        uderiv = chebdifft(uout, 1);
        
        options=optimset('Display','none');
%         options.TolFun = 1e-16;
%         options.TolX = 1e-16;
        % left zero
        zL = fsolve( @(x) chebint(uderiv, x/L), xL, options);
        % right zero
        zR = fsolve( @(x) chebint(uderiv, x/L), xR, options);
        % peak distance (backwards since Chebyshev)
        dist = [dist (zL-zR)];
        
        % difference in peak heights
        ht_diff = chebint(uout, zL) - chebint(uout, zR);
        ht_diffs = [ht_diffs ht_diff];
        
        
    end
end