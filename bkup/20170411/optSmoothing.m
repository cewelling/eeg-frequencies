% Find the optimal smoothing kernel for OscFreq vs. numSwitches correlation

clearvars

global kSize
global pValues
global rValues
%global cutOff

pValues = [];
rValues = [];
for Kernelsize = 0:200:2000
%for cutOff = 0.65:0.05:0.9
    kSize = Kernelsize + 1;
    analysisController
    indexAnal2
    pValues = [pValues P(5,32)];  
    rValues = [rValues R(5,32)];
end

sizes = 0:200:2000;
%cutOffs = 0.65:0.05:0.9;
figure
plot(sizes, pValues, 'ko');
%plot(cutOffs, pValues, 'ko');
xlabel('kernel length (data points)')
%xlabel('filter cut-off (Hz)')
ylabel('OscFreq vs. numSwitches p value')
title('Effect of smoothing on correlation (p-values): Stratus')
%title('Effect of filtering on correlation (p-values): Cumulus')

figure 
%plot(cutOffs, rValues, 'bo')
plot(sizes, rValues, 'bo')
%xlabel('filter cut-off (Hz)')
xlabel('kernel length (data points)')
ylabel('OscFreq vs. numSwitches r value')
title('Effect of smoothing on correlation (r-values) : Stratus')
%title('Effect of filtering on correlation (r-values) : Cumulus')
