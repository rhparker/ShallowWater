%% makes a plot of oscillations near pulse
% u is *half-wave* in this case
% speed c is parameter, which is appended to u
% cutoff is where to cut off plot on the right side

function [decay, freq] = osc_plot(x, u, b, c, L, cutoff)

% find the spatial eigenvalues
% roots of (2/15) nu^4 - b nu^2 + c == 0
nu = roots([(2/15) 0 -b 0 c]);

% oscillations frequency is imag(nu)
% oscillations decay with constant re(nu)
decay = abs(real(nu(1)));
freq  = abs(imag(nu(1)));

% where to start and end plot
l_bound = floor(length(x)/200);
r_bound = floor(length(x) * cutoff/L);
y       = x(l_bound:r_bound);

% scale solution by exp(decay) to recover oscillations
uscaled = u(l_bound:r_bound).*exp(decay*y);
umax    = max(uscaled);

% plot along with sine function of same scaling
figure;
plot(y,uscaled,y,umax*sin(y*freq));
axis([ x(l_bound) x(r_bound) -umax umax]);

legend('rescaled solution','sine function')
title(strcat('scale to see oscillations, speed c =  ',num2str(c)))
 
end
