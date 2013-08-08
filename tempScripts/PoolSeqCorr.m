roi = {'CA1', 'CA3'};
prePost = {'pre', 'post'};
%poolOut = BatchProcess(@SeqCorrs, 'kenji', 'CA1', 'linear', 1, {'prePost'},'pool', 1, struct('poolVar', 'regress'))
colors = {[0,0,0], [1,0,0]};
cntr = 0;
scale = 50;
for kk = 1 : 2
    for m = prePost
        pout = load(['~/data/analysis/kenji/PoolSeqCorr.', roi{kk}, '.', char(m), '.mat']);
        poolOut = pout.out;
        %idx = true(size(poolOut.poolArray, 1), 1);
        idx = poolOut.poolArray(:, 7) <.05 & poolOut.poolArray(:, 7) ~= 0;
        x = [poolOut.poolArray(idx, 2), poolOut.poolArray(idx, 4)];
        pVal = poolOut.poolArray(idx,7);
        cntr = cntr + 1;
        
        if kk == 1
            if strcmp(m, 'pre')
                hSctr(cntr) = scatter(x(:, 1), x(:, 2),'o','filled','SizeData', (1-pVal)*scale + eps, 'CData', colors{kk});
            else
                hSctr(cntr) = scatter(x(:, 1), x(:, 2), '^', 'filled','SizeData', (1-pVal)*scale + eps, 'CData', colors{kk});
            end
        else
            if strcmp(m, 'pre')
                hSctr(cntr) = scatter(x(:, 1), x(:, 2), 'o','filled','SizeData', (1-pVal)*scale + eps, 'CData', colors{kk});
            else
                hSctr(cntr) = scatter(x(:, 1), x(:, 2), '^', 'filled','SizeData', (1-pVal)*scale + eps, 'CData', colors{kk});
            end
        end
        hold on;
    end
end
%ylim([-.015, .015])
set(gca, 'XTick', 1:5);
ylim([ -0.0031    0.0020])
xlim([0, 5.5])
line(xlim, [0, 0], 'Color', 'k');
xlabel('Peak')
ylabel('Slope')
legend(hSctr, {'CA1, Pre', 'CA1, Post', 'CA3, Pre', 'CA3, Post'})



% %%
% be = linspace(min(x(:,2)), max(x(:, 2)), 11);
% [cnt, ~] = histc(x(:, 2), be);
% cnt = cnt ./ sum(cnt);
% h = bar(be(1:end-1), cnt(1:end-1), 'FaceColor', 'k', 'barwidth', .3);
% hold on
% xlabel('Slope')
