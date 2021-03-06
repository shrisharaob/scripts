for kBase = 1 : length(selectedFBs)
    kTrials = matches(strcmp(matches(:, 1), selectedFBs{kBase}), 2);
    nTrs = sum(~cellfun(@isempty, y.cmnCntrs), 2)
    refTr = 1;
   
    for lTr = 1 : nTrs(kBase)
        gt = GenericTrial(selectedFBs{kBase}, kTrials{lTr});
        if lTr == 1, commonClus = gt.LoadCommonClus; end
        kClu = y.clu{kBase}(ismember(y.clu{kBase}, commonClus));
        gt = gt.LoadPF;
        srm = reshape(gt.pfObject.smoothRateMap(:,:, ismember(gt.pfObject.acceptedUnits, kClu)), [], length(kClu));
        tempPk{lTr} = [];
        for mCell = 1 : length(kClu);
            CellIdx = Sub2Ind([length(gt.pfObject.xBin), length(gt.pfObject.yBin)], y.cmnCntrs{kBase, lTr}.cntrPeaks{mCell});
            kSrm = srm(:, mCell);
            mCellPk{mCell} = kSrm(mCellIdx);
            tempPk{lTr}  = [tempPk{lTr}; mCellPk{mCell}];
        end
    end
    pk = cell2mat(tempPk);
    clear tempPk
    for mTr = 1: length(kTrials)
        if mTr ~= refTr
            plot(pk(:, refTr), pk(:, mTr), '*');
            hold all;
        end
    end
end
keyboard;
%end