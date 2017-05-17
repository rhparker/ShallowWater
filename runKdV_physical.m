%% KdV time iteration
% time stepping is Crank-Nicolson / Adams Bashforth 2
% Fourier collocation for space

% load initial data

function [data, time, xnew] = runKdV_physical(x, u, config, iter, sep, options)

    % how often to save data
    if isfield(options, 'save_steps')
    	save_steps = options.save_steps;
    else
        save_steps = 1;
    end

    % default parameters
    par.c = u(end);                         % wave speed
    N_initial = length(x);                  % current grid size
    L_initial = ceil(abs(x(1)));            % current domain length;
    
    N = N_initial;
    L = L_initial;
    h = 2*L/N;  
    
    % current grid spacing
    if strcmp(config.method,'Chebyshev')
        N = N + 2;
        N_initial = N_initial + 2;
        % need third Neumann BC on R
        config.form = 'nonintegrated';
    end
    
    % if we like, we can change stuff here
    
    % if we change stuff, run through Newton solver
    if N ~= N_initial || L ~= L_initial || par.c ~= u(end)
        % interpolate onto a larger grid, or with a longer domain, or with
        % different c
        [xnew, unew] = fsolveequation(x, u, par, N, L, config, 10000);
    else
        % if we don't interpolate
        xnew = x; unew = u;
    end

    uwave = unew(1:end-1);
    
    % compute differentiation matrices
    if strcmp(config.method,'Fourier');
        [D, D2, D3, D4, D5] = D_fourier(N, L);
    elseif strcmp(config.method,'Chebyshev');
        [D, D2, D3, D4, D5] = D_cheb(N, L, config);
    % finite differences
    else
        % grid spacing
        h = 2*L/N;
        [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC, 3);
    end

    %% do timestepping

    % parameters for time-stepping
    k = 0.01;                       % time step size
    total_iter = iter;               % number of iterations 

    uin = uwave;

    % modify wave

    % if sep is an integer, then space out given number of grid points
    if mod(sep, 1) == 0 
        if sep > 0
            um = stretch_wave(uin, sep);
            uin = um;
        elseif sep < 0
            um = compress_wave(uin, abs(sep));
            uin = um;
        end
    % if not an integer, only do for positive sep    
    else
        if sep > 0
            fine_points = 65536;
            fine_spacing = 2*L/fine_points;                      % current grid spacingn n b
            ufine = interpft(uin, fine_points);
            stretch_pts = floor( sep / fine_spacing);
            ustretch = stretch_wave(ufine, stretch_pts);
            um = interpft(ustretch, N);
            uin = um;
        end 
    end

    % if we modify c, wave no longer stationary solution
    % so translates one way or the other based on difference
    % from par.c; this is expected behavior
    c = par.c;
    
    IMEX = true;
    implicit = false;
    
    LN = D5 - D3 + c*D;                 % linear part of spatial operator

    if IMEX
        EB = inv(eye(length(xnew)) - k*LN);            % implicit Euler operator
        CN = inv(eye(length(xnew)) - (k/2)*LN);        % implicit part of C-N operator

        % run one iterator of EB-EF so we can seed CN-AB
        u1 = uin;
        u2 = EB*( uin - k*2*uin.*(D*uin) );
    else
        u1 = uin;
    end

    data = [ uin ];
    time = [ 0 ];

    %% run timestepping scheme
    
%     % data vector for EF/EB
%     dataE = data;

%     wb = waitbar(0,'please wait...');
    
    for iter = [1:total_iter]

        % EF / EB splitting, doesn't appear to work
    %     u = dataE(:,end);
    %     f1 = -2*u.*(D*u);
    %     exp_part = u + k*f1;
    %     unew = EB * exp_part;
    %     dataE = [ dataE unew ];

        if IMEX
        % CN / AB splitting, looks better (hopefully it really is!)
            f1 = -2*u2.*(D*u2);
            f0 = -2*u1.*(D*u1);
            exp_part = u2 + k * ( (3/2)*f1 - (1/2)*f0 + (1/2)*LN*u2 ); 
            unew = CN * exp_part;
            u1 = u2;
            u2 = unew;
        elseif implicit
            uold = data(:,end);
            unew = fsolve( @(u) u - uold - k/2*( KdV(xnew, u,par,config,D,D2,D3,D4,D5,0) + KdV(xnew,uold,par,config,D,D2,D3,D4,D5,0) ), uold);
        end
        
        % if we are saving our data
        if mod(iter, save_steps) == 0
            data = [ data unew ];
            time = [ time k*iter ];
            save output data time;
            disp({iter, total_iter});
%            waitbar(iter / total_iter);
        end
        
    end

%     close(wb);
    
% %   plot final step
%     plot(xnew, uwave, xnew, data(:,1), xnew, data(:,end));
%     legend('original 2-pulse', 'start (initial perturbation)', 'end (after timestepping');
%     title(strcat('KdV5 starting with double pulse, steps: ',num2str(total_iter)));

end

% stretches wave in center by duplicating center value
function uout = stretch_wave(uin, spacing)
    if mod(length(uin), 2) == 0
        center = length(uin)/2 + 1;
    else
        center = (length(uin) + 1)/2;
    end
    center_val = uin(center);
    center_spacer = center_val * ones(spacing + 1, 1);
    left_piece = uin(1 + ceil(spacing/2):center - 1);
    right_piece = flip(left_piece(1:end));
%     if mod(spacing, 2) == 0
%         right_piece = flip(left_piece(2:end));
%     else 
%         right_piece = flip(left_piece(1:end));
%     end
    uout = [left_piece; center_spacer; right_piece ];
    if length(uout) > length(uin)
        uout = uout( 1: length(uin) );
    end
end

function uout = compress_wave(uin, spacing)
    center = length(uin)/2 + 1;
    left_spacer = uin(1) * ones(spacing, 1);
    left_piece = [left_spacer; uin(1:center - spacing)];
    uout = [left_piece; flip(left_piece(2:end-1))];
end



