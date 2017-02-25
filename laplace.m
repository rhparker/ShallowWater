% Laplace equation in 1D with different BCs
% model system to make sure everything is working

% which BCs to use
% periodic = true;

% grids = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384];
grids = [1024, 2048, 4096, 8192, 16384];

nz_eig = zeros([1 length(grids)]);

% which eigenvalue to look at
index = 20;
if periodic
    true_eig = floor(index/2)^2;
else
    true_eig = (index - 1)^2;
end

for i = 1:length(grids) 
    N = grids(i);

    if periodic
        % for periodic BCs, use domain R/2pi, i.e. [-pi, pi]
        % periodic BCs, either Fourier or finite diff
        x = linspace(-pi,pi,N)';  
        x = x(1:end-1);
        h = 2*pi/(N-1);
    
        % % Fourier spectral methods; essentially no error
        % [D, D2] = D_fourier(N-1, pi);
    
        % Finite difference/periodic
        % second order version (to compare with Neumann)
        [D, D2] = D_fdiff(N-1, h, 'periodic', 1);
    else
        % Neumann BCs
        N = grids(i);
        x = linspace(0,pi,N)';  
        h = pi/(N-1);
        [D, D2] = D_fdiff_Neumann(N, h);
    end
    
    % Laplace operator is just second derivative
    % equation is -Laplace
    LN = -D2;
    
%     % run eig on this; requires a full matrix; can be slow
%     [V, LD] = eig(full(LN));
%     lambda = sort(diag(LD));
    
    % run eigs on this, centered at 0.1; faster
    num = 30;
    center = 0.1;
    opts.tol = 10^(-10);
    [V, lambda_D] = eigs(LN, num, center, opts);
    lambda = sort(diag(lambda_D));
    
    % eigenvalue to test
    nz_eig(i) = lambda(index);
end;

% error plots for polynomial order
% use only for fdiff since Fourier converges exp really fast
% so we hit machine accuracy quickly
nz_diff = abs( nz_eig - round(nz_eig(end)) );
nz_diff = abs( nz_eig - true_eig );
xplot = log(grids);
yplot = log(nz_diff);
figure;
hold on;
bestfit = fit(xplot', yplot', 'poly1');
scatter(xplot, yplot, 45, 'filled');
plot(xplot, bestfit(xplot));
title(strcat('Log error vs Log gridsize, order = ',num2str(abs(bestfit.p1))));
xlabel('Log gridsize');
ylabel('Log error');



