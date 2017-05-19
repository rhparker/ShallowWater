% % make nice plots


% plot of first four double pulses
% load 100F;
load 100C;

figure;
plot(xout, ud_out(1:end-1, :));
ylabel('u(x)'); xlabel('x');
legend('2(2)', '2(3)', '2(4)', '2(5)');

figure;
plot(xout, uout(1:end-1));
ylabel('u(x)'); xlabel('x');
