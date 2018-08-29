clear all; close all; clc;
% setFigProps

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pars = [105;106;109;114;115;117;118;119;121;122;123;125;127;128]; % darts rival %190

% pars = [102;103;105;106;109;114;115;117;119;121;123;125;128;130;131;132]; % more strict %190

% pars = [105;106;109;114;115;117;118;119;121;122;123;125;127;128]; % darts sim

pars = [105;106;109;114;115;118;122]; % marz rival (Why did Jackson select these?)

runType = 'rival'; % 'rival' or 'sim'
stimType = 'marz'; % 'darts' or 'marz'

plotWindow = 3; %amount of time to show in plots
epochHalf = 5; %this shouldn't change

numPars = length(pars);

omniTransitions = [];
omniTransitions.f1 = [];
omniTransitions.f2 = [];
omniTransitions(1).type = 'low to high';
omniTransitions(2).type = 'high to low';
omniTransitions(3).type = 'mixed to high';
omniTransitions(4).type = 'mixed to low';
omniTransitions(5).type = 'high to mixed';
omniTransitions(6).type = 'low to mixed';

omniDurations = [];
omniDurations.durs = [];
omniDurations(1).type = 'Slow';
omniDurations(2).type = 'Fast';
omniDurations(3).type = 'Mixed';


% for each subject ID in pars variable, load in subject transition and
% duration data - store in group 'omni variable'
for iPar = 1:length(pars)
    
    par = pars(iPar);
    parName = ['stratus' num2str(par)];
    dataName = ['./epochResults/' parName '_' stimType '_' runType '_Transitions'];
    load(dataName);
    
    for iTrans = 1:length(transitions)
        if size(transitions(iTrans).f1,1) >= 1
            omniTransitions(iTrans).f1 = [omniTransitions(iTrans).f1; nanmean(transitions(iTrans).f1,1)];
            omniTransitions(iTrans).f2 = [omniTransitions(iTrans).f2; nanmean(transitions(iTrans).f2,1)];
        end
    end
    
    omniDurations(1).durs = [omniDurations(1).durs; durations(1).duration];
    omniDurations(2).durs = [omniDurations(2).durs; durations(2).duration];
    omniDurations(3).durs = [omniDurations(3).durs; durations(3).duration];
end


%plotting

%settings
lineClr = ['r';'g'];
lineProps = [];
lineProps.col = cellstr(lineClr);

%loop for transition type

multiplot = figure;
for iCase = 1:size(omniTransitions,2)
    
        %take the mean and error across slow and fast frequency power time
        %courses
        meanF1trace = nanmean(omniTransitions(iCase).f1,1);
        errorF1trace = ste(omniTransitions(iCase).f1);

        meanF2trace = nanmean(omniTransitions(iCase).f2,1);
        errorF2trace = ste(omniTransitions(iCase).f2);
    
        figure(multiplot)
        if length(errorF1trace) > 1 %exclude this plot given insufficient examples
            switch runType
                case 'rival' % if looking at rival runs, assign 6 subplots
                    if iCase < 2 || iCase > 5
                        subplot(2,3,iCase)
                    elseif iCase == 2 || iCase == 3
                        subplot(2,3,(iCase+2))
                    elseif iCase == 4 || iCase == 5
                        subplot(2,3,(iCase-2))
                    end
                case 'sim' % if looking at sim runs, assign 2 subplots
                    if iCase == 1
                        subplot(2,3,2)
                    elseif iCase == 2
                        subplot(2,3,5)
                    end
            end
            title(omniTransitions(iCase).type)
            hold on
            mseb([(-epochHalf+(1/512)):(1/512):epochHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],lineProps,1);
            if iCase < 3 || iCase > 4
                axis([-plotWindow plotWindow -1 1])
            else
%                 axis([(-plotWindow-1) (plotWindow+1) -1 1]) % make mixed
%                 plot wider
                axis([-plotWindow plotWindow -1 1])
            end
            vline(0,'k') %(0,'k','transition')
            xtick = get(gca,'XTick');
            set(gca,'FontSize',10)
            
            %Create single plot
            figure(iCase+1)
            title(omniTransitions(iCase).type)
            hold on
            mseb([(-epochHalf+(1/512)):(1/512):epochHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],lineProps,1);
            if iCase < 3 || iCase > 4
                axis([-plotWindow plotWindow -1 1])
            else
                axis([(-plotWindow-1) (plotWindow+1) -1 1])
            end
            vline(0,'k') %(0,'k','transition')
            xtick = get(gca,'XTick');
            set(gca,'FontSize',10)
            xlabel('Time from Button Press (s)');
            ylabel('Frequency Band Power');
            
            %save single plot
            saveName = ['Stratus_Group_' omniTransitions(iCase).type '-' runType '-' stimType];
            saveas(gcf,['plots/transitions/Group/' saveName],'pdf')
        end
end

figure(multiplot)
suptitle(['Stratus Group Transitions (N = ' num2str(numPars) '): ' runType '-' stimType])
[ax,h1] = suplabel('Time from Button Press (s)','x',[.1 .1 .84 .84]); %[.08 .08 .84 .84]
[ay,h2] = suplabel('Frequency Band Power','y',[.15 .15 .6 .6]);
set(h1,'fontsize',15)
set(h2,'fontsize',15)

saveName = ['Stratus_Group_Transitions-' runType '-' stimType];
saveas(gcf,['plots/transitions/Group/' saveName],'pdf')



%Button press data
slowDur = mean(omniDurations(1).durs);
slowError = std(omniDurations(1).durs);

fastDur = mean(omniDurations(2).durs);
fastError = std(omniDurations(2).durs);

if strcmp(runType,'rival')
    mixedDur = mean(omniDurations(3).durs);
    mixedError = std(omniDurations(3).durs);

    figure
    hold on
    bar(1:3,[slowDur fastDur mixedDur])
    Labels = {'5.67', '8.5 Hz', 'Mixed'};
    set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
    errorbar(1:3,[slowDur fastDur mixedDur],[slowError fastError mixedError],'.')
    title(['Stratus Group Durations (N = ' num2str(numPars) '): ' runType '-' stimType])
    ylabel('seconds')
elseif strcmp(runType,'sim')
    figure
    hold on
    bar(1:2,[slowDur fastDur])
    Labels = {'5.67', '8.5 Hz'};
    set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
    errorbar(1:2,[slowDur fastDur],[slowError fastError],'.')
    title(['Stratus Group Durations (N = ' num2str(numPars) '): ' runType '-' stimType])
    ylabel('seconds')
end