%% KdV time iteration
% time stepping is Crank-Nicolson / Adams Bashforth 2
% Fourier collocation for space

% load initial data

function [data, time, xnew] = runKdV_physical(x, u, config, iter, save_steps)

    % how often to save data
    if ~exist('save_steps','var')
    	save_steps = 1;
    end

    % default parameters
    par.c = u(end);                 % wave speed
    N = length(x);                  % current grid size
    L = -x(1);                      % current domain length;
    h = 2*L/N;                      % current grid spacingn n b
    
    N = 1024;
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
    
    % compute differentiation matrices
    Fourier = strcmp(config.method,'Fourier');
    if Fourier
        [D, D2, D3, D4, D5] = D_fourier(N, L);
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

    um = stretch_wave(uin, 2);
%     um = compress_wave(uin, 1);
    
    % um = uwave + gaussmf(xout, [1 -3.5156]) + gaussmf(xout, [1 3.5156]);
    % um = uwave * 0.90;

    uin = um;

    % if we modify c, wave no longer stationary solution
    % so translates one way or the other based on difference
    % from par.c; this is expected behavior
    c = par.c;
    
    IMEX = true;
    implicit = false;
    
    LN = D5 - D3 + c*D;                 % linear part of spatial operator

    if IMEX
        EB = inv(eye(N) - k*LN);            % implicit Euler operator
        CN = inv(eye(N) - (k/2)*LN);        % implicit part of C-N operator

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

    wb = waitbar(0,'please wait...');
    
    for iter = [1:total_iter]
        % advection eq, for testing
        % unew = RK4( data(:,end), k, @(u) D*u );

        % 5th order KdV equation  / RK4 (appears unstable)
        % unew = RK4( data(:,end), k, @(u) KdV(u,par,N,D,D2,D3,D4,D5,0) );

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
            unew = fsolve( @(u) u - uold - k/2*( KdV(u,par,N,D,D2,D3,D4,D5,0) + KdV(uold,par,N,D,D2,D3,D4,D5,0) ), uold);
        end
        
        % if we are saving our data
        if mod(iter, save_steps) == 0
            data = [ data unew ];
            time = [ time k*iter ];
            waitbar(iter / total_iter);
        end
        
    end

%     close(wb);
    
%     plot(xnew, uwave, xnew, data(:,1), xnew, data(:,end));
%     legend('original 2-pulse', 'start (initial perturbation)', 'end (after timestepping');
%     title(strcat('KdV5 starting with double pulse, steps: ',num2str(total_iter)));

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



