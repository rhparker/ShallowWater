%% plots eigenvalues of linearization about various waves
% loads eigenvalues from file since
% eig takes forever to run

% ranges for plots
% standard ranges
xrange      = 1e-1;
xrange_fine = 1e-3;
yrange      = 1e8;
% range to see distinct eigenvalues near real axis
yrange_zoom = 1e-4;

% single pulse, for comparison purposes
load 2usingle;

figure;
plot(xout,uout);
title(strcat('Single pulse, b = 0.5, c =  ',num2str(par.c)) )

figure;
plot(lambda, '.');
title('Eigenvalues for single pulse');
axis([-xrange_fine xrange_fine -yrange yrange]);

% % double pulse
% load 1udouble1a;
% plot_title = 'double pulse 4, between 1st max and 2nd min';
% 
% figure;
% plot(xout,ud_out(1:end-1));
% title(strcat(plot_title,', c =  ',num2str(par.c)) )
% 
% figure;
% plot(xout,ud_out(1:end-1));
% title(strcat(plot_title,', c =  ',num2str(par.c)) )
% axis([-3 3 -2 0.5]);
% 
% figure;
% plot(lambda, '.');
% title(strcat('eigenvalues, ',plot_title,', c =  ',num2str(par.c)) );
% axis([-xrange_fine xrange_fine -yrange yrange]);


% figure;
% plot(xout,V);
% title(strcat('eigenfunctions, ',plot_title,', c =  ',num2str(par.c)) );