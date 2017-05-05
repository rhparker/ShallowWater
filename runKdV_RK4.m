%% KdV time iteration

% load initial data

function [data, time, xnew] = runKdV_RK4(x, u, config, iter, save_steps)

    % how often to save data
    if ~exist('save','var')
    	save_steps = 1;
    end

    % default parameters
    par.c = u(end);                 % wave speed
    N = length(x);                  % current grid size
    L = -x(1);                      % current domain length
    h = 2*L/N;                      % current grid spacing
    
    % fourier frequencies
    delta_k = 2*pi/(N*h);
    k = [0:delta_k:N/2*delta_k,-(N/2-1)*delta_k:delta_k:-delta_k]';

    % if we change stuff, run through Newton solver
    if N ~= length(x) || L ~= -x(1) || par.c ~= u(end)
        % interpolate onto a larger grid, or with a longer domain, or with
        % different c
        [xnew, unew] = fsolveequation(x, u, par, N, L, config, 10000);
    else
        % if we don't interpolate
        xnew = x; unew = u;
    end

    uwave = unew(1:end-1);

    %% do timestepping

    % parameters for time-stepping
    delta_t = 1 / (N^2);           % time step size
    total_iter = iter;               % number of iterations 

    uin = uwave;

    % modify wave

    um = stretch_wave(uin, 2);
%     um = compress_wave(uin, 1);
    
    % um = uwave + gaussmf(xout, [1 -3.5156]) + gaussmf(xout, [1 3.5156]);
    % um = uwave * 0.90;

    uin = um;

    % if we modify c, wave no longer stationary solution
    % so translates one way or the other based on difference
    % from par.c; this is expected behavior
    c = par.c;

    %% run time stepping
    
    data = [ uin ];
    time = [0];
    
    % FFT of initial data
    F_data = [ fft(uin) ];
    
    uhat = F_data(:,1);
    
    A = 1i * (k.^5 + k.^3 + c*k);
    
    G = -1i*delta_t*k;
    Ehalf = exp(A*delta_t/2); 
    E = Ehalf.^2;
        

    wb = waitbar(0,'please wait...');
    
    for iter = 1:total_iter
   
        % RK4 scheme
        g1 = G.*fft(real( ifft( uhat ) ).^2);
        g2 = G.*fft(real( ifft(Ehalf.*(uhat + g1/2)) ).^2);
        g3 = G.*fft(real( ifft(Ehalf.*uhat + g2/2) ).^2);
        g4 = G.*fft(real( ifft(E.*uhat+Ehalf.*g3) ).^2);
        uhat_next = E.*uhat + (E.*g1 + 2*Ehalf.*(g2+g3) + g4)/6;

        % save only on specified iterations
        if mod(iter, save_steps) == 0
            F_data = [ F_data uhat_next ];
            data = [ data real(ifft(uhat_next)) ];
            time = [ time delta_t*iter ];
            waitbar(iter / total_iter);
        end
        
        uhat = uhat_next;

    end

    close(wb);
    
    plot(xnew, uwave, xnew, data(:,1), xnew, data(:,end));
    legend('original 2-pulse', 'start (initial perturbation)', 'end (after timestepping');
    title(strcat('KdV5 starting with double pulse, steps: ',num2str(total_iter)));

end

% stretches wave in center by duplicating center value
function uout = stretch_wave(uin, spacing)
    center = length(uin)/2 + 1;
    center_val = uin(center);
    center_spacer = center_val * ones(spacing*2 + 1, 1);
    left_piece = uin(1 + spacing:center - 1);
    uout = [left_piece; center_spacer; flip(left_piece(2:end))];
end

function uout = compress_wave(uin, spacing)
    center = length(uin)/2 + 1;
    left_spacer = uin(1) * ones(spacing, 1);
    left_piece = [left_spacer; uin(1:center - spacing)];
    uout = [left_piece; flip(left_piece(2:end-1))];
end



