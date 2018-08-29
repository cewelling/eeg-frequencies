function setElecs( parName )
%setElecs.m Summary of this function goes here
%   Detailed explanation goes here

% load parameters
analysisParams

%% dartboards

% space to store electrodes chosen
elecs = [];

% load files with high SNR electrode information
SNRfiles = dir(['highSNRelecs/' parName '_dart*_' num2str(numElecs) 'elecs.mat']);
if ~isempty(SNRfiles)
    
    elecArray = {}; % store high SNR electrodes from each run
    for thisRun = 1:length(SNRfiles)
        load(['highSNRelecs/' SNRfiles(thisRun).name])
        % only include a run's electrodes if it's highest SNR electrodes exceed the SNR cutoff
%         if nanmean(maxSNRs(2,1:commonElecs)) > minSNR
             elecArray = [elecArray maxSNRs];
%         end
    end
    
    if ~isempty(elecArray)
        % find electrodes with high SNR in all runs
        interElecs = elecArray{1}(1,:);
        for i = 2:length(elecArray)
            interElecs = intersect(interElecs, elecArray{i}(1,:));
        end
        
        % choose the best electrodes among those
        if length(interElecs) > commonElecs
            for i = 1:length(interElecs)
                avgSNR = NaN;
                for j = 1:length(elecArray)
                    %ind = find(ismember(interElecs(i), elecArray{j}(1,:)));
                    avgSNR = nanmean([avgSNR elecArray{j}(2, elecArray{j}(1,:) == interElecs(1,i))]);
                    %avgSNR = nanmean([avgSNR elecArray{j}(2, ind)]);
                end
                interElecs(2,i) = avgSNR; % store average SNRs in row 2
            end
            interElecs(2,isnan(interElecs(2,:))) = -inf;
            [~,sortingIndices] = sort(interElecs(2,:),'descend');
            elecIndices = sortingIndices(1:2);
            elecs = interElecs(1,elecIndices);
        elseif length(interElecs) == commonElecs
            elecs = interElecs;
        else
            error('not enough common electrodes for %s darts', parName)
        end
    end
    save(['highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs'], 'elecs');
end

%% marz

% space to store electrodes chosen
elecs = [];

% load files with high SNR electrode information
SNRfiles = dir(['highSNRelecs/' parName '_marz*_' num2str(numElecs) 'elecs.mat']);
if ~isempty(SNRfiles)
    
    elecArray = {}; % store high SNR electrodes from each run
    for thisRun = 1:length(SNRfiles)
        load(['highSNRelecs/' SNRfiles(thisRun).name])
        % only include a run's electrodes if it's highest SNR electrodes exceed the SNR cutoff
%         if nanmean(maxSNRs(2,commonElecs)) > minSNR
             elecArray = [elecArray maxSNRs];
%         end
    end
    
    if ~isempty(elecArray)
        % find electrodes with high SNR in all runs
        interElecs = elecArray{1}(1,:);
        for i = 2:length(elecArray)
            interElecs = intersect(interElecs, elecArray{i}(1,:));
        end
        
        % choose the best electrodes among those
        if length(interElecs) > commonElecs
            for i = 1:length(interElecs)
                avgSNR = NaN;
                for j = 1:length(elecArray)
                    %ind = find(ismember(interElecs(i), elecArray{j}(1,:)));
                    avgSNR = nanmean([avgSNR elecArray{j}(2, elecArray{j}(1,:) == interElecs(1,i))]);
                    %avgSNR = nanmean([avgSNR elecArray{j}(2, ind)]);
                end
                interElecs(2,i) = avgSNR; % store average SNRs in row 2
            end
            interElecs(2,isnan(interElecs(2,:))) = -inf;
            [~,sortingIndices] = sort(interElecs(2,:),'descend');
            elecIndices = sortingIndices(1:commonElecs);
            elecs = interElecs(1,elecIndices);
        elseif length(interElecs) == commonElecs
            elecs = interElecs;
        else
            error('not enough common electrodes for %s marz', parName)
        end
    end
    save(['highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs'], 'elecs');
end

end

