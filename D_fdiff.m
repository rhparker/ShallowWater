%% Centered finite difference differentiation matrices
% N is number of grid points
% h is mesh size
% BC is boundary conditions ('Neumann', 'periodic')

function [D, D2, D3, D4, D5] = D_fdiff(N, h, BC)

if strcmp(BC,'Neumann')
    [D, D2, D3, D4, D5] = D_fdiff_Neumann(N, h);
else
% we go up to 5th derivative, so use 3-pt stencil
    pts = 3;
    D  = fddiff(N, 1, h, pts, BC);
    D2 = fddiff(N, 2, h, pts, BC);
    D3 = fddiff(N, 3, h, pts, BC);
    D4 = fddiff(N, 4, h, pts, BC);
    D5 = fddiff(N, 5, h, pts, BC);
% otherwise, Neumann BCs (less accurate for now)
end

end