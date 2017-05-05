% compute time derivatives of peak distances

% for now, uses finite differences

function derivs = peak_derivs(t, dist)

    % left and right derivatives, 4th order accuracy
    % we can ignore these if we want
    Lderiv = [-49/20 6 -15/2 20/3 -15/4 6/5 -1/6];
    Rderiv = flip(Lderiv);
    
    % central derivatives
    % 6th order accuracy
%     Cderiv = [1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280];
    halforder = 2;
    weights = Fornberg_weights(0, -halforder:halforder, 1);
    Cderiv = weights(2,:);
   
    N = length(t);
    Dmat = sparse(N, N);
    one_sided_rows = (length(Cderiv) - 1) / 2;
    delta_t = t(2) - t(1);
    
    offset = length(Lderiv) - 1;
    for index = 1:one_sided_rows;
        Dmat(index, index:index + offset) = Lderiv;
        Dmat(end - index + 1, end - index + 1 - offset:end - index + 1) = Rderiv; 
    end

    offset = length(Cderiv) - 1;
    for index = 1 : (N - 2*one_sided_rows)
        Dmat(index+one_sided_rows, index:index+offset) = Cderiv;
    end
    
    Dmat = Dmat / delta_t;
    
    derivs = Dmat * dist;
end