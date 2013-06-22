for ii = 1 : length(cluAB)
    if isempty(cluAB{ii}), cluAB{ii} = []; end 
end
allCntrs{1} = outCntrs{1}(~cellfun(@isempty, cluAB), :);
allCntrs{2} = outCntrs{2}(~cellfun(@isempty, cluAB), :);
cluAB(cellfun(@isempty, cluAB)) = [];

for lArena = 1 : length(allCntrs)
    cntrs = allCntrs{lArena};
    nTrs = sum(~cellfun(@isempty, cntrs), 2);
    for lBase = 1 : size(cntrs, 1)
        for kTr = 1 : nTrs(lBase)
            nCells = length(cntrs{lBase, kTr}.cntrPeaks);
            if nCells > 1
                cellPairs = nchoosek(1 : nCells, 2);
                for mPr = 1 : size(cellPairs, 1) 
                    pkA = cntrs{lBase, kTr}.cntrPeaks{cellPairs(mPr, 1)}; % peaks of cell A
                    pkB = cntrs{lBase, kTr}.cntrPeaks{cellPairs(mPr, 2)}; 
                    if ~isempty(pkA) & ~isempty(pkB)
                        nCntrsA = size(pkA, 1);
                        nCntrsB = size(pkB, 1);
                        cntrPairs = GenPairs(1 : nCntrsA, 1 : nCntrsB);
                        kTrPkDist{kTr, mPr} = vnorm(pkA(cntrPairs(:, 1), :) - pkB(cntrPairs(:, 2), :), 2);
                    else
                        kTrPkDist{kTr, mPr} = nan;
                    end
                end
            else
                kTrPkDist{kTr} = [];
            end
        end
        pkDist{lBase} = kTrPkDist;
        clear kTrPkDist
    end
    keyboard;
    p{lArena} = pkDist;
end