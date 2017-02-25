% takes diff operators and generates exp weight diff operators
% a  :  exponential wt
function [XD, XD2, XD3, XD4, XD5] = D_expwt(D, D2, D3, D4, D5, a)
    % sparse identity matrix; can add to full matrix 
    I = speye(length(D));
       
    XD  = D  - a*I;
    XD2 = D2 - 2*a*D  + (a^2)*I;
    XD3 = D3 - 3*a*D2 + 3*(a^2)*D   - (a^3)*I;
    XD4 = D4 - 4*a*D3 + 6*(a^2)*D2  - 4*(a^3)*D   + (a^4)*I;
    XD5 = D5 - 5*a*D4 + 10*(a^2)*D3 - 10*(a^3)*D2 + 5*(a^4)*D - (a^5)*I;
end