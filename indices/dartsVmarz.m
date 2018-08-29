clearvars

load('dartRival_analMatrix.mat');

domDur = nanmean([analMatrix(:,2), analMatrix(:,3)], 2);
dAnalMatrix = [analMatrix(:, 1) domDur analMatrix(:, 4:6)]
dartsMeans = nanmean(dAnalMatrix);
dartsSte = ste(dAnalMatrix);

load('marzRival_analMatrix.mat');

domDur = nanmean([analMatrix(:,2), analMatrix(:,3)], 2);
mAnalMatrix = [analMatrix(:, 1) domDur analMatrix(:, 4:6)]
marzMeans = nanmean(mAnalMatrix);
marzSte = ste(mAnalMatrix);

% ttests
rejectNull = [];
for i = 1:size(dAnalMatrix, 2)
    h = ttest2(dAnalMatrix(:,i),mAnalMatrix(:,i));
    rejectNull = [rejectNull h];
end

% Bar graphs
figure 
hold on
bar(1:2, [dartsMeans(1), marzMeans(1)])
Labels = {'darts' 'marz'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
set(gca,'fontsize', 11);
errorbar(1:2, [dartsMeans(1), marzMeans(1)], [dartsSte(1) marzSte(1)],'.')
title(['Average Proportion Dominance'])
%ylabel()

figure 
hold on
bar(1:2, [dartsMeans(2), marzMeans(2)])
Labels = {'darts' 'marz'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
set(gca,'fontsize', 11);
errorbar(1:2, [dartsMeans(2), marzMeans(2)], [dartsSte(2) marzSte(2)],'.')
title(['Average Duration of Dominant Percept'])
ylabel('seconds')

figure 
hold on
bar(1:2, [dartsMeans(3), marzMeans(3)])
Labels = {'darts' 'marz'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
set(gca,'fontsize', 11);
errorbar(1:2, [dartsMeans(3), marzMeans(3)], [dartsSte(3) marzSte(3)],'.')
title(['Average Duration of Mixed Percept'])
ylabel('seconds')


