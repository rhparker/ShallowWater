%% Centered finite difference differentiation matrices
% N is number of grid points
% h is mesh size
% BC is boundary conditions ('Neumann', 'periodic')

function [D, D2, D3, D4, D5] = D_fdiff(N, h, BC, stencil)

% 2nd order Neumann BCs, uses older code, i.e. not
%  computed via Fornberg weights
if strcmp(BC, 'Neumann2')
    [D, D2, D3, D4, D5] = D_fdiff_Neumann(N, h);
    
else
    % default is 3-point stencil so we can get to 5th derivative
    if ~exist('stencil','var')
        stencil = 3;
    end

    pts = stencil;
    D  = fddiff(N, 1, h, pts, BC);
    D2 = fddiff(N, 2, h, pts, BC);
    D3 = fddiff(N, 3, h, pts, BC);
    D4 = fddiff(N, 4, h, pts, BC);
    D5 = fddiff(N, 5, h, pts, BC);
end

end