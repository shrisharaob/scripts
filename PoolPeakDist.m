function pkDistAB = PoolPeakDist(varargin)
% pkDistAB = PoolPeakDist(varargin)
% [IF_COMPUTE, datasetType, roi, inArenaPair]
  [IF_COMPUTE, datasetType, roi, inArenaPair] = ...
        DefaultArgs(varargin, {0, 'kenji', 'CA3', {'bigSquare', 'linear'}});

    if IF_COMPUTE
        options.poolVar = 'pkDist';
        if FileExists(['~/data/analysis/kenji/PooledMultiPkDistance', GenFiletag(roi, arena), 'mat'])
            temp = load(['~/data/analysis/kenji/PooledMultiPkDistance', GenFiletag(roi, arena), 'mat'], 'pkDist')
            pkDist = temp.pkDist;
        else
            pkDist = BatchProcess(@MultiPeakPFDistance, 'kenji', roi, inArenaPair, 1, {roi, inArenaPair, [], 1, 0}, 'pool', 0, options);
        end
        fb = pkDist.poolArrayId(:, 1);
        arenas = pkDist.poolArrayId(:, 3);
        filebases = unique(fb);
        pkDistAB = [];
        for lBase = 1 : length(filebases)
            lIdx = cellfun(@strcmp, fb, repmat(filebases(lBase), size(fb)));
            lPoolArrayIdx = pkDist.poolArrayId(lIdx, 4); % row ids of poolarray for lBase
            lArenas = arenas(lIdx);
            arenaAIdx  = cellfun(@strcmp, lArenas, repmat(inArenaPair(1), size(lArenas)));
            arenaBIdx  = cellfun(@strcmp, lArenas, repmat(inArenaPair(2), size(lArenas)));
            poolArrayIdxA = lPoolArrayIdx(arenaAIdx);
            poolArrayIdxB = lPoolArrayIdx(arenaBIdx);
            if sum(arenaAIdx) == 0 | sum(arenaBIdx) == 0, continue; end
            arenaPairs = GenPairs(1 : sum(arenaAIdx), 1 : sum(arenaBIdx));
            for kArenaPr = 1 : size(arenaPairs, 1)
                kArenaPr
                kPkDistA = pkDist.poolArray(cell2mat(poolArrayIdxA(arenaPairs(kArenaPr, 1))), [3]);
                kPkDistB = pkDist.poolArray(cell2mat(poolArrayIdxB(arenaPairs(kArenaPr, 2))), [3]);
                cellPairsA = pkDist.poolArray(cell2mat(poolArrayIdxA(arenaPairs(kArenaPr, 1))), [1, 2]);
                cellPairsB = pkDist.poolArray(cell2mat(poolArrayIdxB(arenaPairs(kArenaPr, 2))), [1, 2]);
                % look for cmn cell pairs
                cmnIdxA = ismember(cellPairsA, cellPairsB, 'rows')
                cmnIdxB = ismember(cellPairsB, cellPairsA, 'rows')
                if sum(cmnIdxA) == 0 | sum(cmnIdxB) == 0, continue; end % if no cmn pairs exist
                cmnCellPairs = GenPairs(find(cmnIdxA), find(cmnIdxB));
                pkDistAB = [pkDistAB; kPkDistA(cmnCellPairs(:, 1)), kPkDistB(cmnCellPairs(:, 2))];
            end
        end
        save(['~/data/analysis/', datasetType, '/', mfilename, GenFiletag(roi, inArenaPair), 'mat'], 'pkDistAB');
    else
        load(['~/data/analysis/', datasetType, '/', mfilename, GenFiletag(roi, inArenaPair), 'mat']);
        figure;
        plot(pkDistAB(:, 1), pkDistAB(:, 2), '*');
        hold on;
        line([min(min(xlim, ylim)), max(max(xlim, ylim))], [min(min(xlim, ylim)), max(max(xlim, ylim))], 'color', 'k');
        axis square;
        xlabel(['pk distance in  ' inArenaPair{1}, ' (pixels)']);
        ylabel(['pk distance in  ' inArenaPair{2}, ' (pixels)']);
        title(roi);
        grid on;
    end

end
           
    