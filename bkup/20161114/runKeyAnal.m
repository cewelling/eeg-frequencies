function [omniData,omniPresent] = runKeyAnal(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,minEpoch)

defaultAnalParams %Load default anal params

if strfind(parName,'caroline') == 1
    date = filenames{1};
    date = date.name;
    [token remain] = strtok(date,'_');
    [date1 remain] = strtok(remain,'_');
    [date2 remain] = strtok(remain,'.');
    [date2 remain] = strtok(date2,'_');
    date = char(datetime(2016,str2num(date1),str2num(date2)));

    if analRunNum < 10
        datafile = [date '_rivalData_caroline0' num2str(analRunNum) '*.txt']; %good one: '14-Mar-2016_rivalData_.txt';
    else
        datafile = [date '_rivalData_caroline' num2str(analRunNum) '*.txt']; %good one: '14-Mar-2016_rivalData_.txt';
    end
    datafile = dir([pressDir datafile]); datafile = datafile.name;
    data = load([pressDir datafile]);
    
elseif strfind(parName,'stratus') == 1
    date = filenames{1};
    date = date.name;
    [token remain] = strtok(date,'_');
    [legend remain] = strtok(remain,'_');
    [runCode remain] = strtok(remain,'_');
    [date1 remain] = strtok(remain,'_');
    [date2 remain] = strtok(remain,'.');
    [date2 remain] = strtok(date2,'_');
    date = char(datetime(2016,str2num(date1),str2num(date2)));
    
    if strfind(runCode,'sim')
        numName = runCode(end);
        numCheck = isstrprop(numName,'digit');
        if numCheck == 0
            numName = [];
        end
        datafile = [date '_simulationData' numName '_' token '.txt']; 
    elseif strfind(runCode,'riv')
        numName = runCode(end);
        datafile = [date '_rivalData_' token '_run' numName '.txt'];
    end
    
    datafile = dir([pressDir datafile]); datafile = datafile.name;
    data = load([pressDir datafile]);
end

for trialIndex = 1:numTrials %1:4 changed from 1:, for 01
    
   trialData = data(find(data(:,1) == trialIndex),:); 
   
   xdata = 1:length(trialData);
   ydata = trialData(:,2);
   timedata = trialData(:,3);
   lastSecond = round(timedata(end));
   timeMultiplier = lastSecond/xdata(end);
   xdata = xdata * timeMultiplier;
   
   figure
   plot(xdata,ydata)
   box off
   presses = ydata > 0;
   try
   axis([0 lastSecond min(ydata(presses)) max(ydata)])
   catch
   end

   
   numPoints = ones(1,length(xdata));
   colorData = repmat([0 0 0]',1,length(xdata));
   
%    figure
%    for index = 1:length(xdata)       
%         if ydata(index) == upArrow
%             scatter(xdata(index),2,'r')
%         elseif ydata(index) == leftArrow
%             scatter(xdata(index),2,'b')
%         elseif ydata(index) == rightArrow
%             scatter(xdata(index),2,'g')
%         else
%             scatter(xdata(index),2,'k')
%         end
%         hold on
%    end
%    axis([0 lastSecond 1.5 2.5])
%    box off
   %xaxis = timedata
    
    
end



%% parse rivalry data
cOmniData = data;
cOmniData = cOmniData(find(cOmniData(:,2) > 0),:); %remove times with no button press
minEventDur = minEpoch; % Have to have pressed the button for (looks like 500 now -Alina) 250 ms to count as a perceptual state

omniStop = []; omniStart = []; omniTime = []; omniKey = [];  omniRow = []; omniTrial = []; omniFreq = [];
for trialIndex = 1:max(cOmniData(:,1)) %changed from 1: for 01

    cData = cOmniData(find(cOmniData(:,1) == trialIndex),:);
    lastSecond = round(cData(end,3));
    
    %make a time vector that is continuous throughout rivalry experiment
    cData(:,6) = cData(:,3) + (trialIndex-1)*lastSecond; 
    
    %start trial on first state press
    cData = cData(min([min(find(cData(:,2)==leftArrow))],[min(find(cData(:,2)==rightArrow))]):end,:);

    elapsedRows = 0;
    for rowIndex = 1:length(cData)-1

        if cData(rowIndex+1,2) == cData(rowIndex,2)
            elapsedRows = elapsedRows + 1;
        else

            keyIndex = cData(rowIndex,2);
            stopTime = cData(rowIndex,3);
            startTime = cData(rowIndex - elapsedRows,3);
            elapsedTime = stopTime - startTime;
            
            %was the person reporting the fast or slow frequency? (
%             freqIndex = 0; %mixed percepts
%             if cData(rowIndex,4) == 1; %fast on left
%                 if keyIndex == leftArrow
%                     freqIndex = 1; %fast frequency
%                 elseif keyIndex == rightArrow
%                     freqIndex = -1; %slow frequency
%                 end
%             elseif cData(rowIndex,4) == 2; %fast on right
%                 if keyIndex == leftArrow
%                     freqIndex = -1; %slow frequency
%                 elseif keyIndex == rightArrow
%                     freqIndex = 1; %fast frequency
%                 end
%             end
            
            freqIndex = 0; %mixed percepts
            if cData(rowIndex,4) == 1 && cData(rowIndex,5) == 1; %fast red on left
                if keyIndex == leftArrow
                    freqIndex = -1; %fast frequency
                elseif keyIndex == rightArrow
                    freqIndex = 1; %slow frequency
                end
            elseif cData(rowIndex,4) == 1 && cData(rowIndex,5) == 2; %slow red on left
                if keyIndex == leftArrow
                    freqIndex = 1; %slow frequency
                elseif keyIndex == rightArrow
                    freqIndex = -1; %fast frequency
                end
            elseif cData(rowIndex,4) == 2 && cData(rowIndex,5) == 1; %fast green on left
                if keyIndex == leftArrow
                    freqIndex = 1; %fast frequency
                elseif keyIndex == rightArrow
                    freqIndex = -1; %slow frequency
                end
            elseif cData(rowIndex,4) == 2 && cData(rowIndex,5) == 2; %slow green on left
                if keyIndex == leftArrow
                    freqIndex = -1; %slow frequency
                elseif keyIndex == rightArrow
                    freqIndex = 1; %fast frequency
                end
            end
            
           
            omniFreq = [omniFreq; freqIndex]; %which frequency was seen
            omniStop = [omniStop; stopTime];  %start time of perceptual episode
            omniStart = [omniStart; startTime]; %end time of perceptual episode
            omniTime = [omniTime; elapsedTime]; %duration of perceptual episode
            omniKey = [omniKey; keyIndex]; %which key was pressed
            omniTrial = [omniTrial; trialIndex]; %which rivalry trial we're on
            omniRow = [omniRow; elapsedRows]; 

            elapsedRows = 0;

            clear stop* start* elapsedTime key*                
        end
    end

end

                    
%final data (identify perceptual states that last for 500ms or longer)
omniData = [omniTrial, omniStart, omniStop, omniTime, omniRow/3, omniKey, omniFreq];
omniData = omniData(find(omniData(:,4) > minEventDur),:)




%save a presentation schedule for the sim runs

if size(cOmniData,2) == 6 %if there is a presentation schedule column
    omniPresent = [];
    
    omniStop = []; omniStart = []; omniTime = []; omniKey = [];  omniRow = []; omniTrial = []; omniFreq = [];
    for trialIndex = 1:max(cOmniData(:,1))
        cData = cOmniData(find(cOmniData(:,1) == trialIndex),:);
        lastSecond = round(cData(end,3));
        
        cData(:,7) = cData(:,3) + (trialIndex-1)*lastSecond;
        
        elapsedRows = 0;
        for rowIndex = 1:length(cData) - 1
            
            if cData(rowIndex+1,6) == cData(rowIndex,6) && rowIndex+1 ~= length(cData)
                elapsedRows = elapsedRows + 1;
            elseif rowIndex+1 == length(cData)
                
                presentCode = cData(rowIndex,6);
                stopTime = cData(rowIndex,3);
                startTime = cData(rowIndex - elapsedRows,3);
                elapsedTime = stopTime - startTime;
                
                if cData(rowIndex,5) == 1 %placement of fast frequency is 'on the left'
                    if presentCode == 1
                        freqIndex = 1; %fast frequency
                    elseif presentCode == 2
                        freqIndex = -1; %slow frequency
                    end
                elseif cData(rowIndex,5) == 2 %place ment of fast frequency is 'on the right'
                    if presentCode == 1
                        freqIndex = -1; %slow frequency
                    elseif presentCode == 2
                        freqIndex = 1; %fast frequency
                    end
                end
                
                omniFreq = [omniFreq; freqIndex]; %which frequency was presented
                omniStop = [omniStop; stopTime];  %start time of presentation
                omniStart = [omniStart; startTime]; %end time of presentation
                omniTime = [omniTime; elapsedTime]; %duration of presentation
                omniTrial = [omniTrial; trialIndex]; %which sim trial we're on
                
                elapsedRows = 0;
                
                clear stop* start* elapsedTime key*
                
            else
                
                presentCode = cData(rowIndex,6);
                stopTime = cData(rowIndex,3);
                startTime = cData(rowIndex - elapsedRows,3);
                elapsedTime = stopTime - startTime;
                
                if cData(rowIndex,5) == 1 %placement of fast frequency is 'on the left'
                    if presentCode == 1
                        freqIndex = 1; %fast frequency
                    elseif presentCode == 2
                        freqIndex = -1; %slow frequency
                    end
                elseif cData(rowIndex,5) == 2 %place ment of fast frequency is 'on the right'
                    if presentCode == 1
                        freqIndex = -1; %slow frequency
                    elseif presentCode == 2
                        freqIndex = 1; %fast frequency
                    end
                end
                
                omniFreq = [omniFreq; freqIndex]; %which frequency was presented
                omniStop = [omniStop; stopTime];  %start time of presentation
                omniStart = [omniStart; startTime]; %end time of presentation
                omniTime = [omniTime; elapsedTime]; %duration of presentation
                omniTrial = [omniTrial; trialIndex]; %which sim trial we're on
                
                elapsedRows = 0;
                
                clear stop* start* elapsedTime key*     
            end
        end
    end
    
    omniPresent = [omniTrial, omniStart, omniStop, omniTime, omniFreq];
    
    
else
    omniPresent = NaN;
end
                    
