function computePLVEEG( iFile, parName, noiseType )
%Calculate Phase Locking Value (PLV)

%Returns the difference ben the phases of the two RLS traces over time
%(where 180 = perfect anti-phase, 0 = perfect in-phase, and 90 = neither)
%and the magnitude of that difference (where 1 = perfect phase locking)
%see demoPLV.m to experiment with sin and cos waves.  To learn more read 
%this blog: https://praneethnamburi.wordpress.com/2011/08/10/plv/ and
%this paper: Lachaux 1999 http://www.ma.utexas.edu/users/davis/reu/ch3/cwt/lachaux.pdf
%CER 2017

% Load parameters
preprocessParams


eegDir = '/Users/robertsonce/Dropbox_MIT/Projects/EEG_Frequencies/analScripts/EEG/rls_data/unsmoothedRLS/';
%cumulus01_dartRival1_principles_Oz

% Generate noise
%generateNoisePLV(sessID,cRunType, noiseType)
load(['indices/PLV/Pilot7b_dartRival_ShuffleNoise.mat'])



%% Load RLS spectra
%maybe this should be smoothed?
if exist([eegDir iFile], 'file')
    load([eegDir iFile])
end

% Initialize variables
maxTrials = size([rls_data.amp],2)/3; % Handle each trial in the RLS data
omniMagPLV = nan(maxTrials,1); 
omniPhasePLV = nan(maxTrials,1);
omniCorr = nan(maxTrials,1);
omniRho = nan(maxTrials,1);
omniP = nan(maxTrials,1);

%% Run analysis
omniRLSlow = [];
omniRLShigh = [];
for iTrial = 1:maxTrials


    %% Subtract low frequency trace from high frequency trace (H-L)
    %Compute difference score of traces

    clear RLSlow RLShigh
    RLSlow = rls_data(2).amp{iTrial};
    RLShigh = rls_data(1).amp{iTrial};
    RLShigh = RLShigh(~isnan(RLShigh));
    RLSlow = RLSlow(~isnan(RLSlow));
    
    t = [rls_time];

    if isempty(RLSlow) || sum(isnan(RLSlow))>100 || sum(isnan(RLSlow)) > 0
        continue
    end

    if isempty(RLShigh) || sum(isnan(RLShigh))>100 || sum(isnan(RLShigh)) > 0
        continue
    end

    %smooth data 
    RLSlow = HRsmoothing(RLSlow, 'gaussian', 501, 71, 0);
    RLShigh = HRsmoothing(RLShigh, 'gaussian', 501, 71, 0);

    minFreq = min(filtSpec.range); %see preprocessingParams
    filtsize = filtSpec.order; %see preprocessingParams


    %% EEG phase-amplitude modulation fitlers

    %detrend
    RLSlow = detrend(RLSlow);
    RLShigh = detrend(RLShigh);

%       %zscore
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
% 
%     derivHi = (sig1b-sig1a)/sampRate;
%     derivLo = (sig2b-sig2a)/sampRate;
% 
%     derivHi = HRsmoothing(derivHi, 'gaussian', 501, 71, 0);
%     derivLo = HRsmoothing(derivLo, 'gaussian', 501, 71, 0);
% 
%     diffDeriv = derivHi - derivLo;
% 
%     [lag_time,twin,xcl]=timewindow_xcorr(derivHi,derivLo,500,0.5,.1,.1,1);

% 
%     midPt = round(size(xcl,2)/2);
%     corVect(:) = xcl(:,midPt);

%         if plotTrialsPLV == 1
%             figure
% 
%             subplot(4,1,1);
%             plot(t,[xlo; xhi]);
%             xlabel(['Phases'])
%             title(['Trial: ' num2str(iTrial) ])
% 
%             subplot(4,1,2);
%             scatter(t,rad2deg(angle(plv))); hold on
%             xlabel(['Diff between phases'])
% 
%             subplot(4,1,3);
%             plot(t,RLSlow); hold on
%             plot(t,RLShigh); hold on
%             xlabel(['Raw Signals'])
% 
%             subplot(4,1,4)
%             plot(t(1:end-1),diffDeriv);
%             xlabel(['Difference of Derivs'])
% 
% 
%         end
%     


    % report stats
    disp(['Trial: ' num2str(iTrial) ' Magnitude of PLV: ' num2str(magPLV)]);
    disp(['Trial: ' num2str(iTrial) ' Phase  of PLV: ' num2str(phasePLV)]);

    omniPhasePLV(iTrial) = [phasePLV];
    omniMagPLV(iTrial) = [magPLV];
    %omniCorr(iTrial) = [nanmean(corVect)];
    
    [rho , p] =corr(RLShigh',RLSlow');
    omniRho(iTrial) = [rho ];
    omniP(iTrial) = [p];
    
    omniRLSlow = [omniRLSlow RLSlow];
    omniRLShigh = [omniRLShigh RLShigh];

end

omniPhasePLV = omniPhasePLV(1:maxTrials);
omniMagPLV = omniMagPLV(1:maxTrials);
omniRho = omniRho(1:maxTrials,:);
omniP = omniP(1:maxTrials,:);

%% Save Data
save(['indices/PLVEEG/' parName '_Signal.mat'], 'omniPhasePLV', 'omniMagPLV','omniPhasePLVNoise','omniMagPLVNoise','omniRho','omniP')

