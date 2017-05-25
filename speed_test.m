% test speed of centers for wave

load 100C;

wave = ud_out_2;
sep = 0.95;

opts.savesteps = 50;
opts.speed = 0;

[data, time] = runKdV_physical(xout, wave, config, 10000, sep, opts);
[dist, hts, ctrs] = peak_dist_cheb(xout, data);
speed = center_speed(time, ctrs)


