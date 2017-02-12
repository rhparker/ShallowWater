%% single wave solution to shallow water equation

function F = shallow_sol(x, b)
    scale = sqrt( 3*(2*b + 1) )/2;
    F = 3*(b + (1/2)) * sech(scale * x).^2;
end