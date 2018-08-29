% practice creating different types of filters and looking at their
% frequency response profile

clearvars 

figure

for cutOff = 0.1:0.1:0.8 %500:500:3000
     Fs = 512; % sampling frequency
    N = 5000; % filter order
    %cutOff = 0.4; % Hz, cutoff frequency
    
    Hd = designfilt('lowpassfir','FilterOrder',N,'CutoffFrequency',cutOff, ...
        'DesignMethod','window','Window','hamming','SampleRate',Fs);

     % For a gaussian kernel:
 %N = 800;
 %X = [-N/2:N/2];
 %norm = normpdf(X, 0, N/2 + 1);

    Npadded = 2^15;
    Y = fft(Hd.Coefficients, Npadded);
    
    %Y = fft(norm, Npadded);
    P2 = abs(Y/Npadded); % take real part
    pow = P2(1:Npadded/2 + 1); % positive frequencies
    pow(2:end-1) = 2*pow(2:end - 1);
    
    freq = Fs*(0:Npadded/2)/Npadded;
    
    %figure
    plot(freq, pow)
    xlim([0 2])
    hold on
    title(['N = ' num2str(N) ' cutOff = 0.1:0.1:0.8 Hz']); %' Cutoff = ' num2str(cutOff) ' Hz']);
end

