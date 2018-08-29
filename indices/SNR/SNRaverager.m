% Average SNRs for each participant and output a convenient list to compare
% with FFT plots

SNRaverages = [];
numParticipants = 1; % initialize
for parNum = 101:137
    parName = ['stratus' num2str(parNum)];
    SNRaverages(numParticipants, 1) = parNum;
    participantSNRs = [];
    for noiseWindow = [6 12]
        for run = 1:3 % 1:3 4:6
            if exist([parName '_run' num2str(run) '_Oz.mat'], 'file')
                load([parName '_run' num2str(run) '_Oz.mat']);
            else
                continue;
            end
            for i = 1:length(meanSNRs)
                % Choose SNRs calculated with the desired frequency and noise window
                if (meanSNRs(i).noiseWindow == noiseWindow)
                    participantSNRs = [participantSNRs; meanSNRs(i).value];
                end
            end
        end
        if noiseWindow == 6
            SNRaverages(numParticipants, 2) = mean(participantSNRs); 
        elseif noiseWindow == 12
            SNRaverages(numParticipants, 3) = mean(participantSNRs);
        end
    end
    numParticipants = numParticipants + 1;
end

csvwrite('SNRaverages_OzDarts.csv', SNRaverages);

              
       
            