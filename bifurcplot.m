load 2double2a_info
load 2double2_info
load 2double1a_info
load 2double1_info

hold on;

% unstable point (saddle) and stable point (center)
scatter(6.9146, 0,'X');
scatter(9.6850, 0, '.');
scatter(12.4514, 0, 'X');
scatter(15.2344, 0, '.');

% plot(smdist1a_1, smderiv1a_1, smdist1a_2, smderiv1a_2, smdist1a_4, smderiv1a_4, smdist1a_6, smderiv1a_6)
% plot(smdist1_5, smderiv1_5, smdist1_4, smderiv1_4,smdist1_3, smderiv1_3);
% plot(dist2_1, deriv2_1, dist2_2, deriv2_2, dist2_4, smderiv2_4, dist2_5, smderiv2_5);
plot(dist2a_1, deriv2a_1, dist2a_2, deriv2a_2, dist2a_4, deriv2a_4, dist2a_6, deriv2a_6);
% plot(smdist1a_14_05, smderiv1a_14_05)
title('phase portrait: derivative of peak distance vs peak distance');
xlabel('peak distance');
ylabel('derivative of peak distance');