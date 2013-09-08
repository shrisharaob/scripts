function PairWiseDist(varargin)

    [datasetType, roi, arena] = DefaultArgs(varargin, {'kenji', 'CA3' ,{'bigSquare', 'linear'}});

    switch datasetType
      case 'kenji'
        %loads outCntrs
            load(['~/data/analysis/', datasetType, '/Contours', GenFiletag(roi, arena), 'mat']);
    end
    for mArena = 1 : length(arena)
        inCntrs = cntrs{mArena}; 
        [tcntrs, selectedClus, allFbs, processedFbs] = StableCntrs(inCntrs(:,1), roi, arena{mArena});
        nBases = size(tcntrs, 1);
        nTrs = sum(~cellfun(@isempty, tcntrs), 2);
        nTrials{mArena} = nTrs;
        clu2select{mArena} = selectedClus;
        outCntrs{mArena} = tcntrs;
        clear tcntrs;
    end
%    for mArena = 1 : length(arena)
      for kBase = 1 : length(allFbs)
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
                            %cluA = cntrs{kBase, 1}
                            if ~isempty(cluA) & ~isempty(cluB)
                                cluAB{kBase} = intersect(cluA, cluB);
                            end
                        end
                    end
                end
            end
      end
        %   end
    for mArena = 1 : length(arena)
        switch datasetType
          case 'kenji'
            
            load(['~/data/analysis/', datasetType, '/Contours', GenFiletag(roi, arena{mArena}), 'mat']);
        end
        inCntrs = cntrs;
        clear cntrs;
        [cntrs, selectedClus, allFbs, processedFbs] = StableCntrs(inCntrs, roi, arena{mArena});
        nBases = size(cntrs, 1);
        nTrs = sum(~cellfun(@isempty, cntrs), 2);
        nTrials{mArena} = nTrs;
        for lBase = 1 : nBases
            for kTr = 1 : nTrs(lBase)
                switch datasetType
                  case 'kenji'
                    load(['~/data/analysis/', datasetType, '/', allFbs{lBase}, '/',  allFbs{lBase}, GenFiletag(roi, arena{mArena}), 'commonClus.mat']);
           %          if length(clu2select{mArena}{lBase}) ~= length(cntrs{lBase, kTr}.cntrPeaks), 
%                         keyboard;
%                         cntrs{lBase, kTr}.cntrPeaks = [];
%                         cluAB{lBase} = [];
%         else
                        cntrs{lBase, kTr}.cntrPeaks(~ismember(clu2select{mArena}{lBase}, cluAB{lBase})) = [];
                        %                   end
                end
                nCells = length(cluAB{lBase});
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
    keyboard
end


