smooth_ksize = 9;
smooth_sd =0.5;

signal_in = HRsmoothing(signal_in, 'gaussian', smooth_ksize, smooth_sd, 0);
