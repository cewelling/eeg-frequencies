clearvars

files = dir('DartsAllElectrodes*');
files = extractfield(files, 'name')';

dartsSNRs = [];
for i = 1:length(files)
    load(files{i});
    dartsSNRs = [dartsSNRs; SNRs'];
end

files = dir('MarzAllElectrodes*');
files = extractfield(files, 'name')';

marzSNRs = [];
for i = 1:length(files)
    load(files{i});
    marzSNRs = [marzSNRs; SNRs'];
end

dartsAverage = (mean(dartsSNRs))';
dartsStd = std(dartsSNRs);
marzAverage = (mean(marzSNRs))';
marzStd = std(marzSNRs);

averages = [dartsAverage marzAverage];
bar(averages);
xlabel('electrodes');
set(gca,'fontsize', 11);
ylabel('SNR: averaged across participants and frequencies');
legend('darts','marz');


