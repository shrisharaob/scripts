function out = SeqTemplate(gt, varargin)
% out = SeqTemplate(gt)
% [IF_COMPUTE, IF_PLOT, smootherSpan]

    [IF_COMPUTE, IF_PLOT, smootherSpan, rateThresh] = DefaultArgs(varargin,{0, 0, 33, 2});
    format short;
    if ~FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat']), IF_COMPUTE = 1; end
    if IF_COMPUTE
        %% linearize position 
        if FileExists([gt.paths.analysis, gt.filebase, '.' gt.trialName, '.linPos.mat'])
            load([gt.paths.analysis, gt.filebase, '.' gt.trialName, '.linPos.mat']);
            linVel = diff(smthLinPos);
        else
            fprintf('\n computing linear pos ...');
            xy = sq(gt.position(:, 1, :));
            [th, linPos] = cart2pol(xy(:, 1), xy(:,2));
            smthLinPos = smooth(linPos, smootherSpan);
            linVel = diff(smthLinPos);
            smthLinPos = smthLinPos(1: end-1);
            save([gt.paths.analysis, gt.filebase, '.' gt.trialName, '.linPos.mat'], 'linPos', 'smthLinPos');
        end
        %% get heading direction'
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
        fwdXY = SelectPeriods(smthLinPos, fwdTimePeriods, 'c');
        rvrsXY = SelectPeriods(smthLinPos, rvrsTimePeriods, 'c');
        [res, clu] = gt.LoadStateRes('RUN', 1, gt.trackingSampleRate, [], 1);
        cluIds  = unique(clu);
        [fwdRes, fwdResIdx] = SelectPeriods(res, fwdTimePeriods, 'd');
        fwdClu = clu(fwdResIdx);
        [rvrsRes, rvrsResIdx] = SelectPeriods(res, rvrsTimePeriods, 'd');
        rvrsClu = clu(rvrsResIdx);
        fwdRM = Compute1DRatemap(fwdRes, fwdClu, fwdXY, gt.trackingSampleRate, 1);
        rvrsRM = Compute1DRatemap(rvrsRes, rvrsClu, rvrsXY, gt.trackingSampleRate, 1);

        %% sort cell id 
        validCellId = cluIds(~cellfun(@isempty, fwdRM));
        fwdRatemaps = cell2mat(fwdRM(find(~cellfun(@isempty, fwdRM))));
        rvrsRatemaps = cell2mat(rvrsRM(find(~cellfun(@isempty, rvrsRM))));
        [maxVal, fwdMaxrm] = max(fwdRatemaps, [], 2);
        validCellId(maxVal < rateThresh) = [];
        fwdMaxrm(maxVal < rateThresh) = [];
        fwdRatemaps(maxVal < rateThresh, :) = [];
        [~, sortedCellIdx] = sort(fwdMaxrm);
        fwdRatemaps = fwdRatemaps(sortedCellIdx, :);
        fwdSortedClu =  validCellId(sortedCellIdx);

        validCellId = cluIds(~cellfun(@isempty, rvrsRM));
        [maxVal, rvrsMaxrm] = max(rvrsRatemaps, [], 2);
        validCellId(maxVal < rateThresh) = [];
        rvrsMaxrm(maxVal < rateThresh) = [];
        rvrsRatemaps(maxVal < rateThresh, :) = [];
        [~, sortedCellIdx] = sort(rvrsMaxrm);
        rvrsSortedClu =  validCellId(sortedCellIdx);
        rvrsRatemaps = rvrsRatemaps(sortedCellIdx, :);
        
        out.fwdRatemaps = fwdRatemaps;
        out.fwdSortedClu = fwdSortedClu;
        out.rvrsRatemaps = rvrsRatemaps;
        out.rvrsSortedClu = rvrsSortedClu;
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat'], 'out');
        out.linPos = linPos;
        out.smthLinPos = smthLinPos;
    else
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat'], 'out');
        load([gt.paths.analysis, gt.filebase, '.' gt.trialName, '.linPos.mat']);
        out.linPos = linPos;
        out.smthLinPos = smthLinPos;
        
    end

    %% display
    if IF_PLOT
        h = figure;
        set(h, 'Position', [1, 25, 1680, 920]);
        subplot(1, 2, 1)
        imagesc(out.fwdRatemaps ./  repmat(max(out.fwdRatemaps, [], 2), 1, size(out.fwdRatemaps, 2))) 
        set(gca, 'ydir', 'normal')   
        set(gca, 'YTick', [1 : length(out.fwdSortedClu)] );
        set(gca, 'YTickLabel', out.fwdSortedClu);
        ylabel('clu id', 'FontSize', 16);
        set(gca, 'FontSize', 12)
        text(repmat(2, size(out.fwdSortedClu)), get(gca, 'Ytick'), num2str(max(out.fwdRatemaps, [], 2), '%.1f'), 'Color', 'w', 'FontWeight', 'bold');
        title('fwd ratemaps');

        subplot(1, 2, 2)
        imagesc(out.rvrsRatemaps ./  repmat(max(out.rvrsRatemaps, [], 2), 1, size(out.rvrsRatemaps, 2)))
        set(gca, 'YTick', [1 : length(out.rvrsSortedClu)] );
        set(gca, 'YTickLabel', out.rvrsSortedClu);
        ylabel('clu id', 'FontSize', 16);
        set(gca, 'FontSize', 12)
        set(gca, 'ydir', 'normal')
        text(repmat(2, size(out.rvrsSortedClu)), get(gca, 'Ytick'), num2str(max(out.rvrsRatemaps, [], 2), '%.1f'), 'Color', 'w', 'FontWeight', 'bold');
        title('reverse ratemaps');

%         % fwd rate maps in reverse cell order
%         subplot(1, 4, 3)
%         [~, idx ] = ismember(out.rvrsSortedClu, out.fwdSortedClu);
%         idx(idx == 0) = [];
%         imagesc(out.fwdRatemaps(idx, :) ./  repmat(max(out.fwdRatemaps(idx, :), [], 2), 1, size(out.fwdRatemaps, 2)));
%         set(gca, 'YTick', [1 : length(idx)] );
%         set(gca, 'YTickLabel', out.fwdSortedClu(idx));
%         set(gca, 'ydir', 'normal')
%         title('fwdRM in rvrs order');
%         text(repmat(2, size(out.fwdSortedClu(idx))), get(gca, 'Ytick'), num2str(max(out.fwdRatemaps(idx, :), [], 2), '%.1f'), 'Color', 'w', 'FontWeight', 'bold');

%         % reverse rm in forward order
%         subplot(1, 4, 4)
%         [~, idx ] = ismember(out.fwdSortedClu, out.rvrsSortedClu);
%         idx(idx == 0) = [];
%         imagesc(out.rvrsRatemaps(idx, :) ./  repmat(max(out.rvrsRatemaps(idx, :), [], 2), 1, size(out.rvrsRatemaps, 2)));
%         set(gca, 'YTick', [1 : length(idx)] );
%         set(gca, 'YTickLabel', out.rvrsSortedClu(idx));
%         set(gca, 'ydir', 'normal')
%         title('rvrsRM in fwd order');
%         text(repmat(2, size(out.rvrsSortedClu(idx))),get(gca, 'Ytick'), num2str(max(out.rvrsRatemaps(idx, :), [], 2), '%.1f'), 'Color', 'w', 'FontWeight', 'bold');

        [~, roi] = SearchKenji(gt.filebase);
        inRegion = [];
        if any(cell2mat(cellfun(@strcmp, roi, repmat({'CA3'}, 1, length(roi)), 'uniformoutput', 0)))
            inRegion = [inRegion, '  CA3  ']; end
        
        if any(cell2mat(cellfun(@strcmp, roi, repmat({'CA1'}, 1, length(roi)), 'uniformoutput', 0)))
            inRegion = [inRegion, '  , CA1']; end
        reportfig(h, mfilename, 0, [gt.filebase, '--', gt.trialName, ' :::: ', inRegion]); 
        close(h); 
        %check if the same cell fires at the same location in forward and reverse dirs
%         figure;
%         r1 = rvrsMaxrm(ismember(out.rvrsSortedClu, out.fwdSortedClu));
%         r2 = rvrsMaxrm(ismember(out.fwdSortedClu, out.rvrsSortedClu));
%         plot(r1, r2, '*');
    end
end


