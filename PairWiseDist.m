function PairWiseDist(varargin)

    [datasetType, roi, arena] = DefaultArgs(varargin, {'kenji', 'CA3' ,{'bigSquare', 'linear'}});

    switch datasetType
      case 'kenji'
        %        if strcmp(arena{1})
        load(['~/data/analysis/kenji/Contours.CA3.bigSquare.mat']); % loads cntrs
        cntrsSquare = cntrs;
        load(['~/data/analysis/kenji/Contours.CA3.linear.mat']);
        cntrsLinear = cntrs;
    end
    
    %[cntrs, selectedClus, allFbs, processedFbs] = StableCntrs(cntrs, roi, arena);
    nBases = size(cntrs, 1);
    nTrs = sum(~cellfun(@isempty, cntrs), 2);
    for lBase = 1 : nBases
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
                    end
                end
            end
        end
        if size(kTrPkDist, 2) == 1
            pkDist{lBase} = cell2mat(kTrPkDist);
        else
            pkDist{lBase} = cell2mat(kTrPkDist');
        end
        clear kTrPkDist
    end
keyboard;

    for kBase = 1 : length(allFbs)
        trPairs = nchoosek(1 : nTrs(kBase), 2);
        for mTrPr = 1 : size(trPairs, 1)
            %            plot(pkDist
        end
    end
end
    
    
