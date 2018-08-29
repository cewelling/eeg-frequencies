% HRsmoothing
% 2012-02-13
function signal_out = HRsmoothing(signal_in, kernal_name, kernal_size, kernal_std, figureSwitch)

% clc
% clear all
% load output1
% signal_in = tosavedata(:,3);
% kernal_name = 'gaussian';
% kernal_size = 7; % currently this parameter must be odd integers >= 1
% kernal_std = 1.2;
% figureSwitch = 1;

kernal_size_repeated = repmat(kernal_size,[1,2]);
signal_out = smoothn(signal_in, kernal_size_repeated, kernal_name, kernal_std);
% signal_out = smoothts(signal_in, kernal_name(1), kernal_size, kernal_std);
if figureSwitch
    signal_length = length(signal_in);
    figure1 = figure();
    freq = 0:(2*pi)/signal_length:pi;
    xdft = fft(signal_in);
    ydft = fft(signal_out);
    subplot(2,1,1); plot(freq,abs(xdft(1:length(signal_in)/2+1)));
    hold on;
    plot(freq,abs(ydft(1:length(signal_in)/2+1)),'r','linewidth',2);
    legend('Original Signal','Bandpass Signal');
    hold off;

    t = ((1:1:signal_length)-3)/4;
    subplot(2,1,2); plot(t,signal_in);
    hold on;
    plot(t,signal_out,'r','linewidth',2);
    legend('Original Signal','Bandpass Signal');
    hold off;
    pause;
    close(figure1);
end