% Find the optimal smoothing kernel for OscFreq vs. numSwitches correlation

clearvars

%global kSize
global pValues
global rValues
global cutOff
global groupNum

%**************************************************************************
groupNum = 0; % 0 for stratus, 1 for cumulus
% Make sure you set the  analysisParams group variable to match groupNum %
%**************************************************************************

pValues = [];
rValues = [];
%for Kernelsize = 0:200:2000
for cutOff = 0.2:0.05:0.8
    %kSize = Kernelsize + 1;
    analysisController
    indexAnal2
    pValues = [pValues P{groupNum + 1}(5,17)];
    rValues = [rValues R{groupNum + 1}(5,17)];
end

%sizes = 0:200:2000;
cutOffs = 0.2:0.05:0.8;
figure
%plot(sizes, pValues, 'ko');
plot(cutOffs, pValues, 'ko');
%xlabel('kernel length (data points)')
xlabel('filter cut-off (Hz)')
ylabel('OscFreq vs. numSwitches p value')

if groupNum + 1 == 1
    title('Effect of smoothing on correlation (p-values): Stratus')
else
    title('Effect of filtering on correlation (p-values): Cumulus')
end

figure
plot(cutOffs, rValues, 'bo')
%plot(sizes, rValues, 'bo')
xlabel('filter cut-off (Hz)')
%xlabel('kernel length (data points)')
ylabel('OscFreq vs. numSwitches r value')

if groupNum + 1 == 1
    title('Effect of smoothing on correlation (r-values) : Stratus')
else
    title('Effect of filtering on correlation (r-values) : Cumulus')
end
