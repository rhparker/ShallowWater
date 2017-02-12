% scale range to [-1, 1]

function scaled_x = unit_scale(x)
    M = max(abs(x));
    scaled_x = x/M;
end