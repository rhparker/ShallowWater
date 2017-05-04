function norms = L2norms(x, data)
    
    norms = [];

    data_size = size(data);
    for index = 1:data_size(2);
        u = data(:,index);
        
        L2norm = trapz(x, u.*u);
        norms = [norms; L2norm];
    end

end