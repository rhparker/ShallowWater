function norms = L2norms(x, data)
    
    norms = [];

    data_size = size(data);
    for index = 1:data_size(2);
        u = data(:,index);
        
        % take abs value in case we have chebyshev points
        % which are written backwards
        L2norm = sqrt( abs( trapz(x, u.*u) ) );
        norms = [norms; L2norm];
    end

end