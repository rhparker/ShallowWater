function dist = peak_distances(x, data)
    dist = [];
    data_size = size(data);
    for index = 1:data_size(2)
        [pks, locs] = findpeaks(data(:, index));
        [d1, d2] = sort(pks);
        pklocs = locs(d2(end-1:end));
        xlocs = x(pklocs);
        pkdist = abs( xlocs(2) - xlocs(1) );
        dist = [dist; pkdist];
    end
end