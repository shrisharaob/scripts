function out = PoolRemapping(datasetType, roi, arena, sparsityThresh, poolVar, IF_PLOT)

if IF_PLOT
out = BatchProcess(@Remapping, datasetType, roi, arena, 0, {datasetType, roi, arena, [],[],[], sparsityThresh}, 'pool', 1, struct('poolVar', poolVar))
save(['~/data/analysis/', datasetType, '.Remapping.', poolVar , GenFiletag(arena, roi) 'mat'], 'out');
figH = figure;
plot(out.poolArray(:,1), out.poolArray(:,2),'.k')
hold  on
maxXY = max(max(xlim, ylim));
xlim([0, maxXY])
ylim([0, maxXY])
switch poolVar
  case 'pkRateMat'
    tag = 'Rate Peak '
    filename =  ['~/thesis/figures/report/RateRemapping_', poolVar, '_', roi, '_', arena];
  case 'pkDistMat'
    tag = 'PF Peak Distances'
    filename =  ['~/thesis/figures/report/', poolVar, '_', roi, '_', arena];
end
xlabel([tag, ' Trial A'])
ylabel([tag, ' Trial B'])
line(xlim, ylim, 'color', 'k')
[r,p, rl, ru] = corrcoef(out.poolArray);
text(6, maxXY - .2 * maxXY, ['r = ', num2str(r(1, 2), '%.2f')], 'FontSize', 12, 'FontWeight', 'b')
%textbp(['r = ', num2str(r(1, 2), '%.2f')], 'FontSize', 12, 'FontWeight', 'b')

ProcessFigure(figH, filename, 6, [2.2, 2.2]);
else
load(['~/data/analysis/', datasetType, '.Remapping.', poolVar , GenFiletag(arena, roi) 'mat'], 'out');
keyboard;
end
keyboard;