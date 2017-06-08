% phase portrait plots

function phaseportrait

    load 100F;

    x = [ 2.05 2.1 2.105 2.2 2.3 5 7 7.5 8 9  ];
    tmax =  [ 220 350 400 100 100 300 800 2500 3000 3000 ];
    
    figure('DefaultAxesFontSize',16);

    hold on;
    
    % plot saddles and center
    c  = uout(end);
    nu = roots([1 0 -1 0 c]);
    a = abs(real(nu(1)));
    b = abs(imag(nu(1)));
    
    scatter(0,  0, 'X');
    scatter(pi/b,  0, '.');
    scatter(2*pi/b, 0, 'X');
    scatter(3*pi/b, 0, '.');
    
    % eigenvalues
    n = [0 1 2 3];
    const = abs(eigenvalues(1,1)) / sqrt(b);
    eigens = const * sqrt(b)*exp(-a * n * pi / (2 * b));

    
    % plot trajectories
    for index = 1:length(x) 
        [T, Y] = find_sol( x(index), tmax(index) );
        plot(Y(:,1), Y(:,2) );
    end

end

function [T, Y] = find_sol(x, tmax)
    options = odeset('RelTol',1e-5,'AbsTol',[1e-5 1e-5]);
    [T,Y] = ode45(@system,[0 tmax],[x 0],options);
end


% ODE system
function dy = system(t,y)
    load 100F;
    
    c  = uout(end);
    nu = roots([1 0 -1 0 c]);
    a = abs(real(nu(1)));
    b = abs(imag(nu(1)));
    
    const = abs(eigenvalues(1,1)) / sqrt(b);
    dy = zeros(2,1);    % a column vector
    dy(1) = y(2);
    dy(2) = const * exp(-a * y(1) ) * sin( b * y(1) );
end