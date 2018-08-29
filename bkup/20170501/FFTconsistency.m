% Test the consistency of FFT frequency analysis data. For data split into
% odd and even runs, do participants correlate with themselves?
%
% Called by:
% Dependencies: groupAnalysisParams.m

% load parameters
groupAnalysisParams

for iGroup = 1:length(groupParNums)
    odds = [];
    evens = [];
    idVect = [];
    for parNum = groupParNums{iGroup}
        numFormat = '%02d';
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        
        % load data
        if exist([indicesDir 'rateFreqs/' parName '_split.mat'],'file')
            load([indicesDir 'rateFreqs/' parName '_split.mat'])
        else
            continue;
        end
        
        idVect = [idVect; parNum];
        
        % collect frequencies calculated from odd trials
        odds = [odds; oscFreq_odd];
        
        % collect frequencies calculated from even trials
        evens = [evens; oscFreq_even];
%         if oscFreq_even > 0.3
%             disp('hi')
%         end
    end
    
    % Do odds and evens correlate?
    [R, P] = corrcoef(odds, evens,'rows', 'pairwise');
    
    % Plot odds against evens
    figure
    plot(odds, evens, 'ko');
    xlabel('odd OscFreqs')
    ylabel('even OscFreqs')
    xax = xlim;
    yax = ylim;
    text(xax(1) + (xax(2) - xax(1))/2, yax(1) + (yax(2) - yax(1))*.95, ['r =' num2str(R(1,2)) '   p =' num2str(P(1,2))])
    
    % plot a trend line
    hold on
    validO = ~isnan(odds);
    validE = ~isnan(evens);
    validBoth = validO & validE;
    p = polyfit(odds(validBoth),evens(validBoth),1);   % p returns 2 coefficients fitting r = a_1 * x + a_2
    line = p(1) .* odds(validBoth) + p(2); % compute a new vector r that has matching datapoints in x
    plot(odds(validBoth), line, '-'); 
    title([analGroupIDs{iGroup} ': Odd/Even OscFreq Consistency'])
end