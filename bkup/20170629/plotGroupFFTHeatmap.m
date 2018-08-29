%%Heatmap of FFTs
%%rows = individuals (sorted by group and rate)
%%columns = frequency
%%intensity = normalized amplitude
%%cer 2017

clear all; close all

%% Set up
%addpath('heatmaps')
%addpath('gramm')

dataDir = 'oscFFTs/';
diags = {'cumulus', 'stratus'};
freqs = load([dataDir 'cumulus01.mat']);
freqs = freqs.freq;
maxFreq = 0.75;
freqs = freqs(1:find(freqs==maxFreq));
maxPars = 30; %max # of participants in either group

%% Process data
omniData = nan(2,maxPars,length(freqs)); %diags x pars x datapoints
normOmniData = nan(2,maxPars,length(freqs)); %diags x pars x datapoints
omniCDF = nan(2,maxPars,length(freqs)); %diags x pars x datapoints
parMax = nan(2,maxPars); %diags x pars 
for diagIndex = 1:2
    
    diagName = diags{diagIndex};
    filelist = dir([dataDir diagName '*.mat']);
    filelist = {filelist.name};
    
    for parIndex = 1:maxPars
        clear parData

        if parIndex > length(filelist)
            continue
        end
        
        parData = load([dataDir filelist{parIndex}]);
        freqs = parData.freq;
        freqs = freqs(1:find(freqs==maxFreq));
        demoData = parData.meanfft;
        cdfData = parData.CDF;
        
        demoData = demoData(1:length(freqs));
        cdfData = cdfData(1:length(freqs));
        
        %normalize and gather
        parMax(diagIndex,parIndex) = (find(demoData == max(demoData)));

        normOmniData(diagIndex,parIndex,1:length(demoData)) = demoData/max(demoData);
        omniData(diagIndex,parIndex,1:length(demoData)) = demoData;
        omniCDF(diagIndex,parIndex,1:length(cdfData)) = cdfData;
       
   
    end
end


%% Plot heatmap figure
[y,sorted] = sort(parMax,2); %sort from low to high within groups
normOmniData = squeeze([normOmniData(1,sorted(1,:),:) normOmniData(2,sorted(2,:),:)]);


parMaxASC = nanmean(parMax(1,:));
parMaxCon = nanmean(parMax(2,:));

for i = 1:size(normOmniData,1)
    if i > maxPars
         rowLabels{i} = 'Controls';
         cind(i) = 1;
    else
         rowLabels{i} = 'Autism';
         cind(i) = 2;
    end
end
keepPars = find(~isnan(normOmniData(:,1)));
normOmniData = normOmniData(keepPars,:);
rowLabels = rowLabels(keepPars);
cind = cind(keepPars);

setFigProps
close all hidden
cmapMatrix = parula(50); % color map
clf
heatmap(normOmniData, freqs,rowLabels); hold on
colorbar(); colormap(cmapMatrix); caxis([0 1])
keepy = ylim;
h = line([parMaxASC,parMaxASC],[keepy(1) keepy(2)/2]);
set(h,'color','r','linewidth',2)
h = line([parMaxCon,parMaxCon],[keepy(2)/2 keepy(2)]);
set(h,'color','r','linewidth',2)

title('Individual Oscillation Frequencies of EEG Data')
box off;
set(gca,'TickDir','out');


%% Plot shaded error figures
omniData = squeeze([omniData(2,:,:) omniData(1,:,:)]);
omniData(find(max(omniData,[],2) > 4),:) = NaN;
omniCDF = squeeze([omniCDF(2,:,:) omniCDF(1,:,:)]);

clear rowLabels
for i = 1:size(omniData,1)
    if i > maxPars
         rowLabels{i} = '1';
         cind(i) = 2;
    else
         rowLabels{i} = '2';
         cind(i) = 1;
    end
end

c=rowLabels;
x=freqs;
y=num2cell(omniData,2);

clear g
g(1,1)=gramm('x',x,'y',y,'color',c);
g(1,2)=gramm('x',x,'y',y,'color',c);

c=rowLabels;
x=freqs;
y=num2cell(omniCDF,2);

g(1,3)=gramm('x',x,'y',y,'color',c);

g(1,1).geom_line();
g(1,1).set_title('Individual FFTs');
%g(1,1).set_color_options('map','matlab'); %It is also possible to provide a custom
% colormap by providing a N-by-3 matrix (columns are R,G,B).

g(1,2).stat_summary();
g(1,2).set_title('Group Averaged FFTs');
%g(1,2).set_color_options('map','matlab');

g(1,3).stat_summary();
g(1,3).set_title('Group Averaged CDFs');
%g(1,3).set_color_options('map','matlab');

g.set_title('Individual Oscillation Frequencies of EEG Data');

figure('Position',[100 100 800 550]);
g.draw();


