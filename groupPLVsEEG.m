%groupPLVs( sessID, cRunType, noiseType )
%Average Phase Locking Value (PLV) over independent sessions
% 
% 
% cRunType = 'dartRival';
%
eegDir = '/Users/robertsonce/Dropbox_MIT/Projects/EEG_Frequencies/analScripts/EEG/rls_data/unsmoothedRLS/';
%cumulus01_dartRival1_principles_Oz

% Load EEG data
dataFiles = dir([eegDir '*stratus*']);
dataFiles = {dataFiles.name};

% Load parameters
preprocessParams

% Set directories
indicesDir = 'indices/';

% Initialize variables
maxSessions = length(dataFiles);
groupMagPLVSignal = [];
groupPhasePLVSignal = [];
groupMagPLVNoise = [];
groupPhasePLVNoise = [];
groupRho=[];
groupP=[];
for sessIndex = 1:maxSessions
    
    iFile = char(dataFiles(sessIndex))
    parName = iFile(1:10);

    computePLVEEG( iFile, parName, 'Shuffle' );
    close all
    
    load([indicesDir 'PLVEEG/' parName '_Signal.mat'])
    
    %if this session doesn't have enough trials left after SNR and motion exclusions
%     if sum(~isnan(omniMagPLV)) < 6
%         continue
%     end
     
    groupMagPLVSignal = [groupMagPLVSignal; omniMagPLV]; %This session's average signal magnitude
    groupPhasePLVSignal = [groupPhasePLVSignal; omniPhasePLV]; %This session's average signal phase
    groupMagPLVNoise = [groupMagPLVNoise; omniMagPLVNoise]; %This session's average noise magnitude
    groupPhasePLVNoise = [groupPhasePLVNoise; omniPhasePLVNoise]; %This session's average noise phase

    groupRho = [groupRho; omniRho]; %
    groupP = [groupP; omniP]; 
    

    clear omni*
    
end
nanmean(groupRho)
nanmean(groupP)

%% Report stats
disp(['Total number of trials across sessions: ' num2str(length(groupMagPLVSignal))]);

groupPhasePLVNoise = groupPhasePLVNoise(~isnan(groupPhasePLVNoise));
groupPhasePLVSignal = groupPhasePLVSignal(~isnan(groupPhasePLVSignal));

muSignal = rad2deg(circ_mean(deg2rad(groupPhasePLVSignal(~isnan(groupPhasePLVSignal))))); %circular mean
rlSignal = rad2deg(circ_r(deg2rad(groupPhasePLVSignal(~isnan(groupPhasePLVSignal))))); %resultant length
muNoise = rad2deg(circ_mean(deg2rad(groupPhasePLVNoise(~isnan(groupPhasePLVNoise))))); %circular mean
rlNoise = rad2deg(circ_r(deg2rad(groupPhasePLVNoise(~isnan(groupPhasePLVNoise))))); %resultant length

disp(['Signal - Circular Mean of PLVs: ' num2str(muSignal)]); %abs OK?
disp(['Signal - Resultant Length of PLVs: ' num2str(rlSignal)]); %abs OK?

disp(['Noise - Circular Mean of PLVs: ' num2str(muNoise)]); %abs OK?
disp(['Noise - Resultant Length of PLVs: ' num2str(rlNoise)]); %abs OK?

%Is signal significantly uniform around 180? Rayleigh's test statistic.
[pUniformitySignal, zUniformitySignal] = circ_rtest(deg2rad(groupPhasePLVSignal(find(~isnan(groupPhasePLVSignal)))));
[pUniformityNoise, zUniformityNoise] = circ_rtest(deg2rad(groupPhasePLVNoise(find(~isnan(groupPhasePLVNoise)))));

%Plot rose histogram
figure
polarhistogram(deg2rad(groupPhasePLVSignal),25,'FaceColor','red','FaceAlpha',.3); hold on;
title(['Signal - Phase Difference of PLV. Rayleigh Test of Uniformity Z = ' num2str(zUniformitySignal) ', p = ' num2str(pUniformitySignal) '  .'])

figure
polarhistogram(deg2rad(groupPhasePLVNoise),25,'FaceColor','blue','FaceAlpha',.3); hold on;
title(['Noise - Phase Difference of PLV. Rayleigh Test of Uniformity Z = ' num2str(zUniformityNoise) ', p = ' num2str(pUniformityNoise) '  .'])

%Significantly greater phase locking magnitude than noise?
[hMag,pMag] = ttest2(groupMagPLVSignal,groupMagPLVNoise)

% [hPhase,pPhase] = ttest2(groupPhasePLV,groupPhasePLVNoise)


figure; 
circ_plot(deg2rad(groupPhasePLVSignal(~isnan(groupPhasePLVSignal))),'hist',[],20,true,true,'linewidth',3,'color','r'); hold on
title('Signal, angles and mean resultant vector')
circ_plot(deg2rad(groupPhasePLVNoise(~isnan(groupPhasePLVNoise))),'hist',[],20,true,true,'linewidth',3,'color','r')
title('Total, angles and mean resultant vector')
legend('Signal','Mean Resultant Vector','Noise')

