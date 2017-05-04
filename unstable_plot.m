load 2double1;

[dataC, xnew] = runKdV(xout, ud_out, config, 3500);

% starting peak distance
start_dist = peak_distances(xnew, [dataC(:,1)]);

[dist, centers] = peak_distances(xnew, dataC, 1);

plot_times = [1, 201, 401, 601, 801]';
plot_x = repmat(xout, 1, length(plot_times));

% plot(plot_x, dataC(:,plot_times));
% legend(cellstr(num2str(plot_times-1)));
% title({'Time stepping for unstable double pulse (times are in plot legend)',
%         strcat('starting distance: ',num2str(start_dist)) });

plot(dist)
title('double pulse peak distance as a function of time step');
xlabel('time step (n)');
ylabel('peak distance');