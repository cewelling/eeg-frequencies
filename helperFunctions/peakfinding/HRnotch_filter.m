function [ y ] = HRnotch_filter( x, q_factor )
%HRnotch_filter remove specific frequency noise
%   Input data x, should be a vector;
%   Output y, same length filtered vector
% 140911 update: 222frames-->search from top 50f, 720frames-->scale to 160
%% parameters
figure_switch = 0;
Fs = 4;
breath_f_range = [0.4,0.6]; % correspond to 24-36 cycles/min
search_range = round(length(x)/222*50); %50; % search from how many most powered frequencies
% q_factor = 5; % 35 default in MATLAB
ab = -3; % -3 default in MATLAB
%% notch filter
% find the breath frequency from top 10 powered frequencies
% it should be within the f range
h = spectrum.periodogram;
msspec = msspectrum(h,x,'Fs',Fs,'NFFT',length(x));
% [y,index] = max(msspec.Data);
% freq = msspec.Frequencies(index);
freq_energy_mat(:,1) =  msspec.Frequencies;
freq_energy_mat(:,2) =  msspec.Data;
freq_energy_mat_st = sortrows(freq_energy_mat,2);
freq_energy_mat_st = flipud(freq_energy_mat_st);
index_within_range = freq_energy_mat_st(1:search_range,1)>breath_f_range(1)&...
    freq_energy_mat_st(1:search_range,1)<breath_f_range(2);
freq = freq_energy_mat_st(find(index_within_range,1,'first'),1);
freq;
% design the filter
% convert that frequency to normalized frequency
W = (2*freq)/Fs;
BW = W/q_factor;
% get the notch filter
[b,a] = iirnotch(W,BW,ab);
% filting
y = filtfilt(b,a,x);
%% show filter results
if figure_switch
    msfilt = msspectrum(h,y,'Fs',Fs,'NFFT',length(x));
    figure1 = figure;
    subplot(221);
    plot(msspec); 
    ylim([-140,-80]);
    title('Before Notch Filtering');
    subplot(223);
    plot(msfilt); 
    ylim([-140,-80]);
    title('After Notch Filtering');

    % figure;
    subplot(222); 
    plot(x);
    xlim([0,222]);ylim([-2,2]*1e-4);
    title('Before Notch Filtering');
    subplot(224);
    plot(y); 
    xlim([0,222]);ylim([-2,2]*1e-4);
    title('After Notch Filtering');
%     save_name = strcat('notch_filter_q',num2str(q_factor,'%02d'),'_ab',...
%         num2str(ab,'%02d'));
%     print(figure1,'-dpng',strcat(save_name,'.png'));
    pause;
    close(figure1);
end
end
