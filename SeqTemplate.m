function out = SeqTemplate(gt)

   %% linearize position 
   fprintf('\n computing linear pos ...');
   xy = sq(gt.position(:, 1, :));
   [th, linearXY] = cart2pol(xy(:, 1), xy(:,2));
   smthLinXY = smooth(linearXY, 64);
   linVel = diff(smthLinXY);
   smthLinXY = smthLinXY(1: end-1);

   %% get heading direction
   fprintf('\n computing heading direction ....');
   minSamples = 220; 
   % forward dir time periods
   idx = InOut(linVel > 0, [], minSamples);
   fwdTimePeriods = idx{1};
   % reverse dir time period
   idx = InOut(linVel < 0, [], minSamples);
   rvrsTimePeriods = idx{1};
   
   %% compute forward and reverse rate maps
   fprintf('\n computing directional rate maps ....');
   fwdXY = SelectPeriods(smthLinXY, fwdTimePeriods, 'c');
   rvrsXY = SelectPeriods(smthLinXY, rvrsTimePeriods, 'c');
   [res, clu] = gt.LoadStateRes('RUN', 1, gt.trackingSampleRate, [], 1);
   cluIds  = unique(clu);
   [fwdRes, fwdResIdx] = SelectPeriods(res, fwdTimePeriods, 'd');
   fwdClu = clu(fwdResIdx);
   [rvrsRes, rvrsResIdx] = SelectPeriods(res, rvrsTimePeriods, 'd');
   rvrsClu = clu(rvrsResIdx);
   fwdRM = Compute1DRatemap(fwdRes, fwdClu, fwdXY, 1);
   rvrsRM = Compute1DRatemap(rvrsRes, rvrsClu, rvrsXY, 1);
   
   %% sort cell id 
   validCellId = cluIds(~cellfun(@isempty, fwdRM));
   fwdRatemaps = cell2mat(fwdRM(find(~cellfun(@isempty, fwdRM))));
   rvrsRatemaps = cell2mat(rvrsRM(find(~cellfun(@isempty, rvrsRM))));
   [~, fwdMaxrm] = max(fwdRatemaps, [], 2);
   [~, sortedCellIdx] = sort(fwdMaxrm);
   fwdRatemaps = fwdRatemaps(sortedCellIdx, :);
   fwdSortedClu =  validCellId(sortedCellIdx);

   validCellId = cluIds(~cellfun(@isempty, rvrsRM));
   [~, rvrsMaxrm] = max(rvrsRatemaps, [], 2);
   [~, sortedCellIdx] = sort(rvrsMaxrm);
   rvrsSortedClu =  validCellId(sortedCellIdx);
   rvrsRatemaps = rvrsRatemaps(sortedCellIdx, :);

   %% display
   subplot(2, 2, 1)
   imagesc(fwdRatemaps ./  repmat(max(fwdRatemaps, [], 2), 1, size(fwdRatemaps, 2))) 
   set(gca, 'YTick', [1 : length(fwdSortedClu)] );
   set(gca, 'YTickLabel', fwdSortedClu);
   ylabel('clu id', 'FontSize', 16);
   set(gca, 'FontSize', 10)
   set(gca, 'ydir', 'normal')   
   
   subplot(2, 2, 2)
   imagesc(rvrsRatemaps ./  repmat(max(rvrsRatemaps, [], 2), 1, size(rvrsRatemaps, 2)))
   set(gca, 'YTick', [1 : length(rvrsSortedClu)] );
   set(gca, 'YTickLabel', rvrsSortedClu);
   ylabel('clu id', 'FontSize', 16);
   set(gca, 'FontSize', 10)
   set(gca, 'ydir', 'normal')

   % fwd rate maps in reverse cell order
   subplot(2, 2, 3)
   [~, idx ] = ismember(rvrsSortedClu, fwdSortedClu);
   idx(idx == 0) = [];
   imagesc(fwdRatemaps(idx, :) ./  repmat(max(fwdRatemaps(idx, :), [], 2), 1, size(fwdRatemaps, 2)));
   set(gca, 'YTick', [1 : length(idx)] );
   set(gca, 'YTickLabel', fwdSortedClu(idx));
   set(gca, 'ydir', 'normal')
   
   % reverse rm in forward order
   subplot(2, 2, 4)
   [~, idx ] = ismember(fwdSortedClu, rvrsSortedClu);
   idx(idx == 0) = [];
   imagesc(rvrsRatemaps(idx, :) ./  repmat(max(rvrsRatemaps(idx, :), [], 2), 1, size(rvrsRatemaps, 2)));
   set(gca, 'YTick', [1 : length(idx)] );
   set(gca, 'YTickLabel', rvrsSortedClu(idx));
   set(gca, 'ydir', 'normal')
   %check if the same cell fires at the same location in forward and reverse dirs
   figure;
   r1 = rvrsMaxrm(ismember(rvrsSortedClu, fwdSortedClu));
   r2 = rvrsMaxrm(ismember(fwdSortedClu, rvrsSortedClu));
   plot(r1, r2, '*');
   keyboard;
end
   
   
