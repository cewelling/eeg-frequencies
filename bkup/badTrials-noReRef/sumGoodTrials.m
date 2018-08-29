% Sum all of the good cumulus trials 
clearvars 

files = dir(['cumulus*']);

trialSum = 0;
for i = 1:length(files)
    load(files(i).name);
    trialSum = trialSum + 18 - sum(sum(excluded));
end

trialSum