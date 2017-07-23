% compute Krein matrix quickly for any z

function [Kz, lambda] = kreinmatrix_z(z)

    % load precalculated data
    load 100F_double2_krein;
    
    nNeg = length(lNeg);
    
    % compute resolvent operator
    Z = M2'*(R2 - z*S2)*M2;
    Z2 = inv( Z(5 + nNeg:end, 5+ nNeg:end) );
    Z(5 + nNeg:end,5 + nNeg:end) = Z2;
    resolv = M2*Z*M2';

    % Cz portion of Krein matrix
    Cz = zeros(nNeg, nNeg);
    for i = 1:nNeg
        for j = 1:nNeg
            Cz(i,j) = (PNegPerp*resolv*PNegPerp*R*vNeg(:,i))' * (PNegPerp*R*vNeg(:,j));
        end
    end
    
    Kz = Rhat - z*ND - Cz;
    lambda = eig(Kz);
end


