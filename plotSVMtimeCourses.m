% Plot shaded error classification timecourses
% Plots timecourses resulting from:
%       % runAirplaneMatrices.m
%       % runSVMtransitions.m
% Dependencies: gramm toolbox

% Choose runID to plot (only edit this)
%%%%%%%%%%%%%%
runID = '005';
%%%%%%%%%%%%%%

accMatDir = ['ML/groupAccMats/'];
timeCourse = dir([accMatDir '*' runID '*']);
load([accMatDir timeCourse.name]);

[token, remain] = strtok(timeCourse.name, '_');
[group, remain] = strtok(remain, '_');

if ~isempty(strfind(remain, 'mixVdom')) || ~isempty(strfind(remain, '3state'))
    withMixed = 1;
else
    withMixed = 0;
end


% % Choose just ONE
% testSubset = 0; % plot for 4 stratuses not used in rest of analysis
% wholeGroup = 1; % plot for whole group
% leftOut = 1;
% 
% withMixed = 1;
% 
% group = 'stratus'; % 'stratus' 'cumulus' 'bothGroups'
% trainSet = 'rivalry'; % rivalry; sim
% testSet = 'rivalTransitions';
% if testSubset
%     runID = '005';
%     latNum = 1;
%     segNum = 1;
% elseif wholeGroup
%     numPars = 15; % Find better way to do this
% elseif leftOut
%     numPars = 33;
% end
% 
% % load the timecourse of interest
% %accMatDir = ['ML/accMats/'];
% if testSubset
%     accMatDir = ['ML/accMats/run' runID '/'];
%     load([accMatDir group '_train-' trainSet '_test-' testSet '_lat' num2str(latNum) '_seg' num2str(segNum)]);
% elseif wholeGroup
%     accMatDir = ['ML/groupAccMats/'];
%     if withMixed
%         load([accMatDir group '_train-' trainSet '_test-' testSet '_' num2str(numPars) 'pars_3stateMixed_numPointsClassThreshold.mat']);
%     else
%         load([accMatDir group '_train-' trainSet '_test-' testSet '_' num2str(numPars) 'pars']);
%     end
% elseif leftOut
%     accMatDir = ['ML/groupAccMats/'];
%     if withMixed
%         load([accMatDir group '_train-' trainSet '_test-' testSet '_' num2str(numPars) 'left_out_pars_withMixed']);
%     else
%         load([accMatDir group '_train-' trainSet '_test-' testSet '_' num2str(numPars) 'left_out_pars_PROP']);
%     end
% end

% Maybe just input the filename here?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(group, 'stratus')
    groupLabel = 'Controls';
elseif strcmp(group, 'cumulus')
    groupLabel = 'ASC';
elseif strcmp(group, 'bothGroups')
    groupLabel = 'All participants';
end

%% Plot accuracies

x=teTime;
%y=num2cell(accuracyMat,2);
y=accuracyMat;

clear g
g(1,1)=gramm('x',x,'y',y);
g(1,2)=gramm('x',x,'y',y);

g(1,1).geom_line();
g(1,1).set_title('Individual Classification Accuracies');
%g(1,1).set_color_options('map','matlab'); %It is also possible to provide a custom
% colormap by providing a N-by-3 matrix (columns are R,G,B).

g(1,2).stat_summary('type','sem','setylim',false);
g(1,2).set_title('Group Averaged Classification Accuracies');

if withMixed
    x = teTime;
    y = domAccMat;
    g(1,3)=gramm('x',x,'y',y);
    %g(1,3).stat_summary('type','sem','setylim',false);
    g(1,3).geom_line();
    g(1,3).set_title('Dominant Classification Accuracies');

    x = teTime;
    y = mixAccMat;
    g(1,4)=gramm('x',x,'y',y);
    %g(1,4).stat_summary('type','sem','setylim',false);
   g(1,4).geom_line();
    g(1,4).set_title('Mixed Classification Accuracies');
end

%g.set_title([groupLabel ': train on sim segments, test on rivalry transitions']);
g.set_title(groupLabel);

% Set appropriate names for axes
g.set_names('x','time from reported transition (s)','y','accuracy');

% Font sizes
g.set_text_options('base_size', 14, 'label_scaling', 1, 'legend_scaling', 1, 'legend_title_scaling', 1, 'title_scaling', 1, 'big_title_scaling', 1.4);

figure('Position',[100 100 800 550]);
g.draw();

%% Plot posterior probabilities

x=teTime;
%y=num2cell(accuracyMat,2);
y=postProbMat;

clear g
g(1,1)=gramm('x',x,'y',y);
g(1,2)=gramm('x',x,'y',y);

g(1,1).geom_line();
g(1,1).set_title('Individual Posterior Probabilities');
%g(1,1).set_color_options('map','matlab'); %It is also possible to provide a custom
% colormap by providing a N-by-3 matrix (columns are R,G,B).

%g(1,2).stat_summary('type','sem','setylim',false);
g(1,2).stat_summary('type','sem','setylim',false);

g(1,2).set_title('Group Averaged Posterior Probabilities');

g.set_title([groupLabel ': train on sim segments, test on rivalry transitions'], 'FontSize', 16);

% Set appropriate names for legends
g.set_names('x','time from reported transition (s)','y','posterior probability');

figure('Position',[100 100 800 550]);
g.draw();

