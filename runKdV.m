%% KdV time iteration

% load initial data

function [dataC, xnew] = runKdV(iter)

    load 2double1a;
    uout = ud_out;

    % default parameters
    par.c = uout(end);              % wave speed
    N = length(xout);               % current grid size
    L = -xout(1);                   % current domain length;
    h = 2*L/N;                      % current grid spacingn n b

    N = 512;
    
    % if we change stuff, run through Newton solver
    if N ~= length(xout) || L ~= -xout(1) || par.c ~= uout(end)
        % interpolate onto a larger grid, or with a longer domain, or with
        % different c
        [xnew, unew] = fsolveequation(xout, uout, par, N, L, config, 10000);
    else
        % if we don't interpolate
        xnew = xout; unew = uout;
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
    k = 0.1;                       % time step size
    total_iter = iter;               % number of iterations 

    uin = uwave;

    % modify wave

    um = stretch_wave(uin, 3);
%     um = compress_wave(uin, 1);
    
    % um = uwave + gaussmf(xout, [1 -3.5156]) + gaussmf(xout, [1 3.5156]);
    % um = uwave * 0.90;

    uin = um;

    % if we modify c, wave no longer stationary solution
    % so translates one way or the other based on difference
    % from par.c; this is expected behavior
    c = par.c;

    LN = D5 - D3 + c*D;                 % linear part of spatial operator
    EB = inv(eye(N) - k*LN);            % implicit Euler operator
    CN = inv(eye(N) - (k/2)*LN);        % implicit part of C-N operator

    % run one iterator of EB-EF so we can seed CN-AB
    unew = EB*( uin - k*2*uin.*(D*uin) );
    data = [ uin unew ];


    %%
    % data vectors for EF/EB and CN/AB
    dataE = data;
    dataC = data;

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

        % CN / AB splitting, looks better (hopefully it really is!)
        u = dataC(:,end-1:end);
        f1 = -2*u(:,2).*(D*u(:,2));
        f0 = -2*u(:,1).*(D*u(:,1));
        exp_part = u(:,2) + k * ( (3/2)*f1 - (1/2)*f0 + (1/2)*LN*u(:,2) ); 
        unew = CN * exp_part;
        dataC = [ dataC unew ];
    end


    plot(xnew, uwave, xnew, dataC(:,1), xnew, dataC(:,end));
    legend('original 2-pulse', 'start (initial perturbation)', 'end (after timestepping');

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



