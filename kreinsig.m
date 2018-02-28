% find Krein signatures

% load the spectrums we want
load 1500F_domain_105;

% parameters

% for c = 150;
target = 3.9858;
cutoff = 0.01;
ess_cutoff = 1;
total = length(output_L);

output_eigens = [];
output_ess = [];
count = length(output_lambda(:,1));
eigen_V = [];
ess_V = [];
kreinfns = [];
kreinsigs = [];

figure('DefaultAxesFontSize',16);
set(gca,'LineWidth',2);
xlabel('Domain length (L)');
ylabel('Eigenvalues');
MarkerSize = 40;
hold on;

for index = 1:total
    x = output_x(:, index);
    u = output_u(:, index);
    par.c = u(end);
    % find eigenvalues in unweighted spectrum
    % only want the one with positive real part
    lambda = output_lambda(:,index);
    eigen_index = find( abs( imag(lambda) - target) < cutoff);
    eigens = lambda(eigen_index);
    output_eigens = [output_eigens eigens];    
%     eigen_V = [eigen_V output_V(:, eigen_index, index)];
    eigen_V = output_V(:, eigen_index, index);

    % find essential spectrum point nearest our eigenvalue
    ess_index = find( abs( imag(lambda) - imag(eigens) )  < ess_cutoff);
    % this list includes our eigenvalues, so we need to remove it
    ess_index   = setdiff(ess_index, eigen_index);
    nearest_ess = lambda( ess_index );
    output_ess = [output_ess nearest_ess];
%     ess_V = [ess_V output_V(:, ess_index, index)];
    ess_V = output_V(:, ess_index, index);

    % plot these values for a check
    plot( output_L(index)*[1], imag(eigens), '.', 'Color', 'red', 'MarkerSize', MarkerSize);
    plot( output_L(index)*[1], imag(nearest_ess), '.', 'Color', 'blue', 'MarkerSize', MarkerSize);
    legend('Interaction Eigenvalue','Essential Spectrum Eigenvalue');
    
    N = length(x);
    L = ceil(abs(x(1)));
    % no symmetry for Jacobian
    config.symmetry = 'none';
    H = get_jacobian(x, u(1:end-1), par, config, 'integrated', 0);
    
    krein = [ (1/2)*eigen_V'*H*eigen_V ; (1/2)*ess_V'*H*ess_V ];
    kreinfns  = [kreinfns krein];
    kreinsigs = [kreinsigs sign(real(krein))];
    
end

