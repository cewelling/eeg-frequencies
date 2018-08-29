


datafile = '21-Mar-2016_rivalData_caroline01_day6.txt'; %good one: '14-Mar-2016_rivalData_.txt';
dataDir = '../../runScripts/Results/';
data = load([dataDir datafile]);

upArrow = 38;
leftArrow = 37;
rightArrow = 39;

for trialIndex = 1:4 %1:4
    
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
   axis([0 lastSecond min(ydata) max(ydata)])

   
   numPoints = ones(1,length(xdata));
   colorData = repmat([0 0 0]',1,length(xdata));
   
   figure
   for index = 1:length(xdata)       
        if ydata(index) == upArrow
            scatter(xdata(index),2,'r')
        elseif ydata(index) == leftArrow
            scatter(xdata(index),2,'b')
        elseif ydata(index) == rightArrow
            scatter(xdata(index),2,'g')
        else
            scatter(xdata(index),2,'k')
        end
        hold on
   end
   axis([0 lastSecond 1.5 2.5])
   box off
   %xaxis = timedata
    
    
end



%% parse rivalry data
cOmniData = data;
cOmniData = cOmniData(find(cOmniData(:,2) > 0),:); %remove times with no button press
minEventDur = 0.25; % Have to have pressed the button for 250 ms to count as a perceptual state

omniStop = []; omniStart = []; omniTime = []; omniKey = [];  omniRow = []; omniTrial = []; omniFreq = [];
for trialIndex = 1:max(cOmniData(:,1))

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
            freqIndex = 0; %mixed percepts
            if cData(rowIndex,4) == 1; %fast on left
                if keyIndex == leftArrow
                    freqIndex = 1; %fast frequency
                elseif keyIndex == rightArrow
                    freqIndex = -1; %slow frequency
                end
            elseif cData(rowIndex,4) == 2; %fast on right
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

                    
%final data
omniData = [omniTrial, omniStart, omniStop, omniTime, omniRow/3, omniKey, omniFreq];
omniData = omniData(find(omniData(:,4) > minEventDur),:)