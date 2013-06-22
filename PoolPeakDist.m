function out = PoolPeakDist(pkDist, varargin)

    [datasetType, roi, inArenaPair] = DefaultArgs(varargin, {'kenji', 'CA3', {'bigSquare', 'linear'}});

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
        arenaPairs = GenPairs(1 : sum(arenaAIdx), 1 : sum(arenaBIdx));
        for kArenaPr = 1 : size(arenaPairs, 1)
            kPkDistA = pkDist.poolArray(cell2mat(poolArrayIdxA(arenaPairs(kArenaPr, 1))), [3]);
            kPkDistB = pkDist.poolArray(cell2mat(poolArrayIdxB(arenaPairs(kArenaPr, 2))), [3]);
            cellPairsA = pkDist.poolArray(cell2mat(poolArrayIdxA(arenaPairs(kArenaPr, 1))), [1, 2]);
            cellPairsB = pkDist.poolArray(cell2mat(poolArrayIdxB(arenaPairs(kArenaPr, 2))), [1, 2]);
            % look for cmn cell pairs
            cmnIdxA = ismember(cellPairsA, cellPairsB, 'rows');
            cmnIdxB = ismember(cellPairsB, cellPairsA, 'rows');
            if sum(cmnIdxA) == 0 | sum(cmnIdxB) == 0, continue; end % if no cmn pairs exist
            cmnCellPairs = GenPairs(find(cmnIdxA), find(cmnIdxB));
            pkDistAB = [pkDistAB; kPkDistA(cmnCellPairs(:, 1)), kPkDistB(cmnCellPairs(:, 2))];
        end
    end
keyboard;
end
           
    