%% Look at time indices defined by button presses
function [lowFreqReportData, highFreqReportData] = epochOfInt(cat,buttonReport, dataStruct, latency,conservative,halfCentroidTime)

    %function takes input variables:
    %cat: string, 'yes' or 'no', to indicate whether to handle concatenated trials or independent trials
    %buttonReport: expects omniData output from runKeyAnal, provides button
    %pressing time points do determine periods of time to consider in EEG
    %data
    %dataStruct: preprocessed EEG data
    %latency: time(s) to shift the reported percept start and end points
    
switch cat

    case 'yes'
        
    concatOmniData = buttonReport;
    
    for i = 1:size(concatOmniData,1)
        if concatOmniData(i,1) > 1
            concatOmniData(i,2) = concatOmniData(i,2) + (30 * (concatOmniData(i,1)-1));
            concatOmniData(i,3) = concatOmniData(i,3) + (30 * (concatOmniData(i,1)-1));
        end
    end
    
    % Time shift
    concatOmniData(:,2) = concatOmniData(:,2) - latency;
    concatOmniData(:,3) = concatOmniData(:,3) - latency;
    
    lowFreqReportData = dataStruct;            
    meanTimes = mean([concatOmniData(:,2) concatOmniData(:,3)],2);
    
    for i = 1:length(lowFreqReportData.time{1,1})
        [distance closeEpochRow] = min(abs(meanTimes-lowFreqReportData.time{1,1}(1,i)));
        if lowFreqReportData.time{1,1}(1,i) > concatOmniData(closeEpochRow,2) && lowFreqReportData.time{1,1}(1,i) < concatOmniData(closeEpochRow,3) && concatOmniData(closeEpochRow,7) == -1
            lowFreqReportData.trial{1,1}(:,i) = lowFreqReportData.trial{1,1}(:,i);
        else
            lowFreqReportData.trial{1,1}(:,i) = 0;
        end
    end
    
    highFreqReportData = dataStruct;            
    meanTimes = mean([concatOmniData(:,2) concatOmniData(:,3)],2);
    
    for i = 1:length(highFreqReportData.time{1,1})
        [distance closeEpochRow] = min(abs(meanTimes-highFreqReportData.time{1,1}(1,i)));
        if highFreqReportData.time{1,1}(1,i) > concatOmniData(closeEpochRow,2) && highFreqReportData.time{1,1}(1,i) < concatOmniData(closeEpochRow,3) && concatOmniData(closeEpochRow,7) == 1
            highFreqReportData.trial{1,1}(:,i) = highFreqReportData.trial{1,1}(:,i);
        else
            highFreqReportData.trial{1,1}(:,i) = 0;
        end
    end
    
    case 'no'
        
    % Time shift
    buttonReport(:,2) = buttonReport(:,2) - latency;
    buttonReport(:,3) = buttonReport(:,3) - latency;
        
    lowFreqReportData = dataStruct;
    highFreqReportData = dataStruct;
    
    meanTimes = mean([buttonReport(:,2) buttonReport(:,3)],2);
    meanTimes = [buttonReport(:,1), meanTimes];
    
    switch conservative
        
        case 'no' %include the edges of the perceptual report
    
            for trial = 1:size(dataStruct.trial,2)

                thisButtonReport = buttonReport(buttonReport(:,1) == trial,:);
                %thisButtonReport = buttonReport(buttonReport(:,1) == (trial+1),:); %stratus13

                for i = 1:length(lowFreqReportData.time{1,trial})
                    [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == trial,2)-lowFreqReportData.time{1,trial}(1,i)));
                    if lowFreqReportData.time{1,trial}(1,i) > thisButtonReport(closeEpochRow,2) && lowFreqReportData.time{1,1}(1,i) < thisButtonReport(closeEpochRow,3) && thisButtonReport(closeEpochRow,7) == -1
                        lowFreqReportData.trial{1,trial}(:,i) = lowFreqReportData.trial{1,trial}(:,i);
                    else
                        lowFreqReportData.trial{1,trial}(:,i) = 0;
                    end
                end


                for i = 1:length(highFreqReportData.time{1,trial})
                    [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == trial,2)-highFreqReportData.time{1,trial}(1,i)));
                    if highFreqReportData.time{1,trial}(1,i) > thisButtonReport(closeEpochRow,2) && highFreqReportData.time{1,trial}(1,i) < thisButtonReport(closeEpochRow,3) && thisButtonReport(closeEpochRow,7) == 1
                        highFreqReportData.trial{1,trial}(:,i) = highFreqReportData.trial{1,trial}(:,i);
                    else
                        highFreqReportData.trial{1,trial}(:,i) = 0;
                    end
                end

            end        
        
        case 'yes' %remove the edges of the report
            
            for trial = 1:size(dataStruct.trial,2)

                thisButtonReport = buttonReport(buttonReport(:,1) == trial,:);
                %thisButtonReport = buttonReport(buttonReport(:,1) == (trial+1),:); %stratus13

                for i = 1:length(lowFreqReportData.time{1,trial})
                    [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == trial,2)-lowFreqReportData.time{1,trial}(1,i)));
%                     if lowFreqReportData.time{1,trial}(1,i) > (distance-(halfCentroidTime*1.5)) && lowFreqReportData.time{1,1}(1,i) < (distance+(halfCentroidTime*0.5)) && thisButtonReport(closeEpochRow,7) == -1
                    if lowFreqReportData.time{1,trial}(1,i) > thisButtonReport(closeEpochRow,2) && lowFreqReportData.time{1,1}(1,i) < (thisButtonReport(closeEpochRow,2)+(0.6*thisButtonReport(closeEpochRow,4))) && thisButtonReport(closeEpochRow,7) == -1
                        lowFreqReportData.trial{1,trial}(:,i) = lowFreqReportData.trial{1,trial}(:,i);
                    else
                        lowFreqReportData.trial{1,trial}(:,i) = 0;
                    end
                end


                for i = 1:length(highFreqReportData.time{1,trial})
                    [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == trial,2)-highFreqReportData.time{1,trial}(1,i)));
%                     if highFreqReportData.time{1,trial}(1,i) > (distance-(halfCentroidTime*1.5)) && highFreqReportData.time{1,trial}(1,i) < (distance-(halfCentroidTime*0.5)) && thisButtonReport(closeEpochRow,7) == 1
                    if highFreqReportData.time{1,trial}(1,i) > thisButtonReport(closeEpochRow,2) && highFreqReportData.time{1,trial}(1,i) < (thisButtonReport(closeEpochRow,2)+(0.6*thisButtonReport(closeEpochRow,4))) && thisButtonReport(closeEpochRow,7) == 1
                        highFreqReportData.trial{1,trial}(:,i) = highFreqReportData.trial{1,trial}(:,i);
                    else
                        highFreqReportData.trial{1,trial}(:,i) = 0;
                    end
                end

            end 
            
    end
    
end