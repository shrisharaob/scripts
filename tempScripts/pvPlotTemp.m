function out = pvPlotTemp(gt, roi, arena)

% [pv, av, dp] = PopVecTimeCourse(gt);
% dp(isnan(dp)) = 0;
% subplot(3, 1, 2)
% plot(dp)
% xlabel('Theta Cycles');
% ylabel('Dot product');


%[pv, av, out.dp] = PopVecTimeCourse(gt, 0, 1, 3);
temp = load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), 'CHUNKS.3.PopVecTimeCourse.mat']);
% hold on;
% plot(dp, '*-', 'Color',)
% ylim([0, 1])
out.dp = temp.out.dotProd;

end



% out = BatchProcess(@pvPlotTemp, 'kenji', 'CA1', 'bigSquare', 1, {}, 'pool', 0, struct('poolVar', 'dp'))
% plot(out.poolArray', 'o-', 'Color',  [0.5,0.5,0.5], 'MarkerSize', 1, 'LineWidth', 0.09)
% hold on;
% plot(mean(out.poolArray), 'x-', 'Color', 'k', 'LineWidth', .6, 'MarkerSize', 2)
% plot(mean(out.poolArray), 's-', 'Color', 'k', 'LineWidth', .6, 'MarkerSize', 2)

% xlim([.9, 3.1])
% ylim([0, 1])
% set(gca, 'xtick', 1:3)
% xlabel('Segment #');
% ylabel('Dot Product');
% %