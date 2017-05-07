function [dist, centers] = peak_distances(x, data, skip)

    if ~exist('scale','var')
    	skip = 1;
    end

    dist = [];
    centers = [];
    data_size = size(data);
    
    % make the fine grid
    finepoints = 20000;
    L = -x(1);
    xfine = linspace(-L, L, finepoints+1);
    xfine = xfine(1:end-1);

    data_size_scaled = floor(data_size(2) / skip);
    
    for index = (1:data_size_scaled)*skip

        % cubic spline method
        uspline = spline(x, data(:, index), xfine); 
        xout = xfine;
        uout = uspline;
        
%         % Fourier interpolation method
%         ufourierinterp = interpft(data(:, index), finepoints);
%         xout = xfine;
%         uout = ufourierinterp;
        
        [pks, locs] = findpeaks(uout);
        [d1, d2] = sort(pks);
        pklocs = locs(d2(end-1:end));
        xlocs = xout(pklocs);
        pkdist = abs( xlocs(2) - xlocs(1) );
        
        center = (xlocs(2) + xlocs(1))/2;
        dist = [dist; pkdist];
        centers = [centers; center];
    end
end