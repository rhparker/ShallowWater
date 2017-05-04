%% make oscillation plots about nonlinear center

load 2double1a;

% number of iterations
k = 0.01;
iter = 1000;

x = xout_1024;
u = ud_out_1024;

[dataC, xnew] = runKdV(x, u, config, iter);

% find peak distance for double pulse
start_dist = peak_distances(xnew, [u(1:end-1)]);

% find peak distances in timestepping 
[dist, centers] = peak_distances(xnew, dataC, 1);

% plot on same plot
time_axis = (1:length(dist))*start_dist;
center_line = start_dist*ones(length(dist), 1);

figure;
plot(time_axis, dist, time_axis, center_line);
title(strcat('Double Pulse peak distance vs time, starting distance: ',num2str(dist(1))));
legend('time stepping','stationary double pulse');
xlabel('time');
ylabel('peak distance');


