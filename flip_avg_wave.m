% generate flipped wave and average of flipped 
% and nonflipped wave

function [uflip, uavg] = flip_avg_wave(u, config)
    
    % generate flipped wave u(-x)
    % recall that u(end) is the wave speed
    c = u(end);
    if strcmp(config.BC, 'periodic')
        uflip  = [u(1) ; flip(u(2:end-1)) ; u(end) ];
    else
       uflip  = [ flip(u(1:end-1)) ; u(end) ];
    end;

    % average original wave and flipped wave
    uavg      = 0.5*(u + uflip);
    uavg(end) = c;
end