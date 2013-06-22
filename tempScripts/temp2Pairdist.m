c1 = cntrs{1};
c2 = cntrs{2};
idx1 = logical(sum(~cellfun(@isempty, c1), 2));
idx2 = logical(sum(~cellfun(@isempty, c2), 2));
idx = idx1 & idx2;
c1(~idx, :) = [];
c2(~idx, :) = [];
nBases = sum(idx);
nTr1 = sum(~cellfun(@isempty, c1), 2);
nTr2 = sum(~cellfun(@isempty, c2), 2);
for kBase = 1 : nBases 
   trPairs = GenPairs(1:nTr1(kBase), 1:nTr2(kBase));
    for lTrPr = 1 : size(trPairs, 1)
        clu1 = c1{kBase, trPairs(lTrPr, 1)}.cluId'
        clu2 = c2{kBase, trPairs(lTrPr, 2)}.cluId'
        clu12 = intersect(clu1, clu2);
        if ~isempty(clu12)
            if length(clu12) > 1
                clu12pairs = nchoosek(1 : length(clu12), 2); % pairs that are active in both 
                for mCellPr = 1 : size(clu12pairs, 1) 
                %% pair dist in 1
                    pkA = c1{kBase, trPairs(lTrPr, 1)}.cntrPeaks{clu12pairs(mCellPr, 1)}; % peaks of cell A
                    pkB = c1{kBase, trPairs(lTrPr, 1)}.cntrPeaks{clu12pairs(mCellPr, 2)}; % peaks of cell A
                    if ~isempty(pkA) & ~isempty(pkB)
                        nCntrsA = size(pkA, 1);
                        nCntrsB = size(pkB, 1);
                        cntrPairs = GenPairs(1 : nCntrsA, 1 : nCntrsB);
                        kTrPkDist{lTrPr, mCellPr} = vnorm(pkA(cntrPairs(:, 1), :) - pkB(cntrPairs(:, 2), :), 2);
                    else
                        kTrPkDist{kTr, mCellPr} = nan;
                    end
                    pkdist1{kBase} = kTrPkDist;
                    %% pair dist in 2
                    pkA = c2{kBase, trPairs(lTrPr, 2)}.cntrPeaks{clu12pairs(mCellPr, 1)}; % peaks of cell A
                    pkB = c2{kBase, trPairs(lTrPr, 2)}.cntrPeaks{clu12pairs(mCellPr, 2)}; % peaks of cell A
                    if ~isempty(pkA) & ~isempty(pkB)
                        nCntrsA = size(pkA, 1);
                        nCntrsB = size(pkB, 1);
                        cntrPairs = GenPairs(1 : nCntrsA, 1 : nCntrsB);
                        kTrPkDist{lTrPr, mCellPr} = vnorm(pkA(cntrPairs(:, 1), :) - pkB(cntrPairs(:, 2), :), 2);
                    else
                        kTrPkDist{kTr, mCellPr} = nan;
                    end
                    pkdist2{kBase} = kTrPkDist;
                end
            end
         
        end
    end
end
keyboard;