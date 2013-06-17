function PairWiseDist(varargin)

    [datasetType, roi, arena] = DefaultArgs(varargin, {'kenji', 'CA3' ,{'bigSquare', 'linear'}});

    for mArena = 1 : length(arena)
        switch datasetType
          case 'kenji'
            
            load(['~/data/analysis/', datasetType, '/Contours', GenFiletag(roi, arena{mArena}), 'mat']);
        end
        inCntrs = cntrs;
        [cntrs, selectedClus, allFbs, processedFbs] = StableCntrs(inCntrs, roi, arena{mArena});
        nBases = size(cntrs, 1);
        nTrs = sum(~cellfun(@isempty, cntrs), 2);
        nTrials{mArena} = nTrs;
        for lBase = 1 : nBases
            for kTr = 1 : nTrs(lBase)
                switch datasetType
                  case 'kenji'
                    load(['~/data/analysis/', datasetType, '/', allFbs{lBase}, GenFiletag(roi, arena{mArena}), 'commonClus.mat']);
                end
                
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
%             if size(kTrPkDist, 2) == 1
%                 pkDist{lBase} = cell2mat(kTrPkDist)';
%             else
%                 pkDist{lBase} = cell2mat(kTrPkDist');
%             end
            pkDist{lBase} = kTrPkDist;
            clear kTrPkDist
        end
 
        %       pairWiseDist{mArena} = pkDist;
       pairWiseDist{mArena} = pkDist;
        clu2select{mArena} = selectedClus;
    end

    %% 
    for kBase = 1 : length(allFbs)
        if kBase == 6; keyboard; end
        nArenas = length(arena);
        if nArenas > 1
            arenaPairs = nchoosek(1 : nArenas, 2);
            for mArenaPr = 1 : size(arenaPairs, 1)
                nTrsA = nTrials{arenaPairs(mArenaPr, 1)}(kBase);
                nTrsB = nTrials{arenaPairs(mArenaPr, 2)}(kBase);
                if nTrsA > 0 & nTrsB > 0
                    abTrialPairs = GenPairs(1 : nTrsA, 1 : nTrsB);
                    for mTrPr = 1 : size(abTrialPairs, 1)
                        cluA = clu2select{arenaPairs(mArenaPr, 1)}{kBase};
                        cluB = clu2select{arenaPairs(mArenaPr, 2)}{kBase};
                        if ~isempty(cluA) & ~isempty(cluB)
                            cluAB{kBase} = intersect(cluA, cluB);
                        end
                    end
                end
            end
        end
    end
end


