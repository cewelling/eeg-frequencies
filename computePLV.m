function computePLV(parName, cRunType, noiseType, paramsFlag)
%Calculate Phase Locking Value (PLV)
%transitions = aggTransitions(parName, iPar, runName, date, EEGfile, paramsFlag);
%Returns the difference ben the phases of the two RLS traces over time
%(where 180 = perfect anti-phase, 0 = perfect in-phase, and 90 = neither)
%and the magnitude of that difference (where 1 = perfect phase locking)
%see demoPLV.m to experiment with sin and cos waves.  To learn more read 
%this blog: https://praneethnamburi.wordpress.com/2011/08/10/plv/ and
%this paper: Lachaux 1999 http://www.ma.utexas.edu/users/davis/reu/ch3/cwt/lachaux.pdf
%CER 2017



%% Set-up
% load parameters
analysisParams

% ensure that SNRs have been calculated for this run
thisRunSNRs = dir(['pre-processing/highSNRelecs/' parName '_' cRunType '*' electrodeSet.name '.mat']);

% load and store SNRs
for iTrial = 1:numTrials
    if exist(['pre-processing/highSNRelecs/' parName '_' cRunType '_trial' num2str(iTrial) '_' electrodeSet.name '.mat'],'file')
        load(['pre-processing/highSNRelecs/' parName '_' cRunType '_trial' num2str(iTrial) '_' electrodeSet.name '.mat']);
        snrVect(iTrial) = maxSNRs(2, ismember(maxSNRs(1,:), snrElecs));
        snrVect_lo(iTrial) = lFreqSNRs(2, ismember(lFreqSNRs(1,:), snrElecs));
        snrVect_hi(iTrial) = hFreqSNRs(2, ismember(hFreqSNRs(1,:), snrElecs));
    else % trial doesn't exist (due to mistake)
        snrVect(iTrial) = NaN;
        snrVect_lo(iTrial) = NaN;
        snrVect_hi(iTrial) = NaN;
    end
end

% Generate noise
% Temporarily just loading noise from mouse data
%generateNoisePLV(sessID,cRunType, noiseType)
load(['/Users/robertsonce/Dropbox/Projects/mouseBinoc/rivalryScripts/analScripts/ephys/indices/PLV/Pilot7b_dartRival_ShuffleNoise.mat'])

% Initialize variables
maxTrials = numTrials*(length(runIndices)); 
omniMagPLV = nan(maxTrials,1); 
omniPhasePLV = nan(maxTrials,1);
omniCorr = nan(maxTrials,1);
omniRho = nan(maxTrials,1);
omniP = nan(maxTrials,1);
omniRLSlow = [];
omniRLShigh = [];
normRLS_data = [];


%% Run analysis
for runIndex = runIndices
    runName = [cRunType num2str(runIndex)];

    % Load RLS spectra
    if exist(['rls_data/smoothedRLS/' parName '_' runName '_' analFreqLabel '_' electrodeSet.name '.mat'], 'file')
        load(['rls_data/smoothedRLS/' parName '_' runName '_' analFreqLabel '_' electrodeSet.name '.mat'])
    else
        
        continue
    end

    for iTrial = 1:size(rls_data(1).amp, 2)

        clear RLSlow RLShigh

        RLSlow = rls_data(2).amp{iTrial};
        RLShigh = rls_data(1).amp{iTrial};
        RLShigh = RLShigh(~isnan(RLShigh));
        RLSlow = RLSlow(~isnan(RLSlow));
        rls_time = rls_time(end-length(~isnan(RLShigh))+1:end);

        t = [rls_time]; 

        % account for trials that weren't recorded properly (mistakes)
        if isempty(RLSlow) || sum(isnan(RLSlow))>100 ||sum((RLSlow)==0)>1000 
            continue
           
        end

        if isempty(RLShigh) || sum(isnan(RLShigh))>100 ||sum((RLShigh)==0)>1000 
            continue
        end

        %smooth data 
        RLSlow = HRsmoothing(RLSlow, 'gaussian', 501, 71, 0);
        RLShigh = HRsmoothing(RLShigh, 'gaussian', 501, 71, 0);

        %% EEG phase-amplitude modulation fitlers
        %detrend
        RLSlow = detrend(RLSlow);
        RLShigh = detrend(RLShigh);

        %zscore
        RLSlow = zscore(RLSlow);
        RLShigh = zscore(RLShigh);
        filteredHi = RLShigh;
        filteredLo = RLSlow;

        %hilbert
        hilbertHi=angle(hilbert(filteredHi));
        hilbertLo=angle(hilbert(filteredLo));

        %rename
        xhi = hilbertHi;
        xlo = hilbertLo;

        %absolute phase difference
        clear i
        pd = xhi-xlo; 

        %put into circular space
        plv=exp(i*pd);

        %mean phase-locking value across the trial
        meanPLV=mean(plv);

        %magnitude of plv
        magPLV = abs(meanPLV);

        %phase of plv
        phasePLV = (rad2deg(angle(meanPLV))); 

        sig1 = filteredHi;
        sig2 = filteredLo;

        %for derivs
        sig1a = sig1(1:end-1);
        sig2a = sig2(1:end-1);
        sig1b = sig1(2:end);
        sig2b = sig2(2:end);
        sampRate;

    %         if plotTrialsPLV == 1
    
%                 figure
%     
%                 subplot(3,1,1);
%                 plot(t,[xlo; xhi]);
%                 xlabel(['Phases'])
%                 title([runName ' Trial: ' num2str(iTrial) ])
%     
%                 subplot(3,1,2);
%                 scatter(t,rad2deg(angle(plv))); hold on
%                 xlabel(['Diff between phases'])
%     
%                 subplot(3,1,3);
%                 plot(t,RLSlow); hold on
%                 plot(t,RLShigh); hold on
%                 xlabel(['Raw Signals'])

    % 
    %         end


        omniPhasePLV(iTrial*runIndex) = [phasePLV];
        omniMagPLV(iTrial*runIndex) = [magPLV];

        [rho,p] =corr(RLShigh',RLSlow');
        omniRho(iTrial*runIndex) = [rho ];
        omniP(iTrial*runIndex) = [p];

        omniRLSlow = [omniRLSlow RLSlow];
        omniRLShigh = [omniRLShigh RLShigh];

    end
end
omniPhasePLV = omniPhasePLV(1:maxTrials);
omniMagPLV = omniMagPLV(1:maxTrials);
omniRho = omniRho(1:maxTrials,:);
omniP = omniP(1:maxTrials,:);

%% Save Data
save(['indices/PLV/'  parName '_' cRunType '_' analFreqLabel '_' electrodeSet.name '_Signal.mat'], 'omniPhasePLV', 'omniMagPLV','omniPhasePLVNoise','omniMagPLVNoise','omniRho','omniP')

