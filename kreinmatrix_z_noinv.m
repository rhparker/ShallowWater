% compute Krein matrix quickly for any z
% this version does not invert the resolvent matrix

function [Kz, lambda] = kreinmatrix_z_noinv(z, filename)

    % load precalculated data
    if exist('filename','var')
        load(filename);
    else
        load 100F_double2_krein;
    end
        
    nNeg = length(lNeg);
    
    % Z is the inverse of the resolvent operator
    % this is the inverse in M2 coordinates
    Z = M2'*(R2 - z*S2)*M2;
    
    % we want to solve Zy = PNegPerp*R*vNeg(:,i))
    % except we only want to do that in M2 coordinates
    
    % Cz portion of Krein matrix
    % since this thing is diagonal, let's just force
    % that to be the case and not even calculate the
    % off-diagonal elements
    Cz = zeros(nNeg, nNeg);
    for i = 1:nNeg
%         for j = 1:nNeg
            y = M2'* PNegPerp*R*vNeg(:,i);
            % solve Zx = y, but only for the (potentially) invertible part
            Z2 = Z(5 + nNeg : end, 5 + nNeg : end);
            y2 = y(5 + nNeg : end);
            x2 = linsolve(Z2, y2);
            x = [zeros(4 + nNeg, 1) ; x2];
            % convert this back to regular coordinates
            x = M2*x;
%             Cz(i,j) = (x' * PNegPerp*R*vNeg(:,j));
            Cz(i,i) = (x' * PNegPerp*R*vNeg(:,i));
%         end
    end
     
%     Kz = Rhat - z*ND - Cz;
    
    % enforce diagonal
    Kz = diag(diag(Rhat)) - z*ND - Cz;
    lambda = diag(Kz);
end


