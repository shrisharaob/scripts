function out = CatTrialCCG(filebase, roi,  arena, state, varargin)
% thi script concatinates trials in the same arena and computes ccgs for overlapping subfields,
% since the subfields might be slightly shifted across trials, in di
    
    [cellPairs, type, datasetType, binSize, maxTimeLag, IF_PLOT] = ...
        DefaultArgs(varargin, {[], 'load', 'kenji', 10e-3, 1000e-3, 1});
    out = [];
    if ~FileExists(['~/data/analysis/', datasetType, '/', filebase, '/' filebase, '.', mfilename, '.mat']), type = 'compute'; end
    switch datasetType
      case 'kenji'
        searchStruct.roi = roi;
        searchStruct.arena = arena;
        matches = SearchKenji(searchStruct);
        filebases = unique(matches(:, 1));
    end
    switch type
      case 'compute'
        [matches, fbRoi] =  SearchKenji(filebase);
        if ~any(cellfun(@strcmp, fbRoi, repmat({roi}, size(fbRoi)))); return; end
        trialNames = matches(cellfun(@strcmp, matches(:, 3), repmat({arena}, size(matches(:, 3)))), 2);
        trialPeriods = [];
        load(['~/data/analysis/', datasetType, '/Contours', GenFiletag(roi, arena), 'mat'], 'cntrs', 'filebases');
        [stableCntrs, cellIds] = StableCntrs(cntrs(cellfun(@strcmp, filebases, repmat({filebase}, size(filebases))), :), roi, arena);
        cellIds = cell2mat(cellIds);
        poolTrCntr = 0;
        for mTr = 1 : length(trialNames)
            gt = GenericTrial(filebase, trialNames{mTr});
            mpf = stableCntrs{mTr};
            if length(trialNames) == 1, cellIds = mpf.cluId; end
            if length(cellIds) > 1, cellPairs = nchoosek(cellIds, 2); else, fprintf('\n no pairs found'); return; end
            cntrVertices = mpf.cntrVertices;
            cntrPeaks = mpf.cntrPeaks;
            trialPeriods = [trialPeriods; gt.LoadStatePeriods(state, 0, 0)];
            [res, clu] = gt.LoadStateRes(state, 1, [], [], 1); % res in xy fs
            clu = clu(ismember(clu, cellIds));
            res = res(ismember(clu, cellIds));
            binnedPos = BinPos(gt);
            nCellPairs = size(cellPairs, 1);
            stsp = gt.LoadStatePeriods(state, gt.lfpSampleRate, 1);
            stsp = ConvertFs(stsp - gt.trialPeriods(1), gt.lfpSampleRate, gt.trackingSampleRate);
            counterPk = 0;
            PF_OVERLAP = false(nCellPairs, 1);
            ccgPars = [];
            ncntrpr(mTr) = 0;
            catRes{1} = [];
            catClu{1} = [];
            cntrprCounter  = 0;
            for mCellPair = 1 : nCellPairs
                cntrA = cntrVertices{cellIds == cellPairs(mCellPair, 1)}; 
                cntrB = cntrVertices{cellIds == cellPairs(mCellPair, 2)};
                nCntrA = length(cntrA);
                nCntrB = length(cntrB);
                if nCntrA >= 1 && nCntrB >= 1, 
                    cntrPairs = nchoosek([1 : nCntrA, 1 : nCntrB], 2); % all pairs of selected sub contours
                    cntrPairs(cntrPairs(:, 1) > nCntrA, :) = [];
                    cntrPairs(cntrPairs(:, 2) > nCntrB, :) = [];
                    cntrPairs = sortrows(unique(cntrPairs, 'rows'));
                    pkA = cntrPeaks{cellIds == cellPairs(mCellPair, 1)};
                    pkB = cntrPeaks{cellIds == cellPairs(mCellPair, 2)};
                    validCntrCnt = 0;
                    ncntrpr(mTr) = ncntrpr(mTr) + size(cntrPairs, 1);
                    tcatRes = [];
                    tcatClu = [];
                    for kCntrPr = 1 : size(cntrPairs, 1)
                        % if the sub cntrs overlap
                        %pkDistAB(validCntrCnt + 1) = norm( pkA(cntrPairs(kCntrPr, 1), :) - pkB(cntrPairs(kCntrPr, 2), :));
                        if any(InPolygon(cntrA{cntrPairs(kCntrPr, 1)}, cntrB{cntrPairs(kCntrPr, 2)}));
                            POS_IN_CNTR_A = InPolygon(binnedPos, cntrA{cntrPairs(kCntrPr, 1)});
                            POS_IN_CNTR_B = InPolygon(binnedPos, cntrB{cntrPairs(kCntrPr, 2)});
                            POS_IN_AB = POS_IN_CNTR_A | POS_IN_CNTR_B;
                            inout = InOut(POS_IN_AB);
                            cntrprCounter = cntrprCounter + 1;
                            timesInAB = ConvertFs(IntersectRanges(inout{1}, stsp), gt.trackingSampleRate, gt.sampleRate);
                            res1 = SelectPeriods(res(clu == cellPairs(mCellPair, 1)), timesInAB, 'd', 1, 1);
                            res2 = SelectPeriods(res(clu == cellPairs(mCellPair, 2)), timesInAB, 'd', 1, 1);
                            mClu = [ones(length(res1), 1) * cellPairs(mCellPair, 1); ones(length(res2), 1) * cellPairs(mCellPair, 2)];
                            if mTr == 1, catRes{cntrprCounter} = []; catClu{cntrprCounter} = []; end
                            if isempty(catRes{cntrprCounter}), catOffset = 0;
                            else, catOffset = catRes{cntrprCounter}(end); end
                                catRes{cntrprCounter} = [catRes{cntrprCounter}; res1 + catOffset ; res2 + catOffset];
                                catClu{cntrprCounter} = [catClu{cntrprCounter}; mClu];
                                catPairs(cntrprCounter, :) = cellPairs(mCellPair, :);
                                pkDistAB(mCellPair, cntrprCounter) = norm( pkA(cntrPairs(kCntrPr, 1), :) - pkB(cntrPairs(kCntrPr, 2), :));
                                poolTrCntr = poolTrCntr + 1;
                                %offvsdist(poolTrCntr, : ) = 
                                PF_OVERLAP(mCellPair) = true;
                        end
                    end
                end
            end
        end
        if exist('pkDistAB', 'var')
            p = pkDistAB(pkDistAB ~= 0);
            dvsoff = [];
            sts = CatTrials(gt.filebase, roi, arena, state, gt.sampleRate);
            %% ccg
            validCntrPrCounter = 0;
            for mCntrPr = 1 : length(catRes)
                if length(unique(catClu{mCntrPr})) > 1
                    validCntrPrCounter = validCntrPrCounter + 1;
                    out.cellPairs(validCntrPrCounter, :) = catPairs(mCntrPr, :);
                    tccg = CCGPars(catRes{mCntrPr}, catClu{mCntrPr}, gt.sampleRate, [], [], binSize, maxTimeLag,  [], 0, IF_PLOT);
                    ccgPars{validCntrPrCounter} = tccg{1};
                else
                    continue;
                end
                %% DISPLAY 
                if IF_PLOT
                    subplot(2, 2, 1); hold on;
                    plot(ccgPars{validCntrPrCounter}.smoothTAx, ccgPars{validCntrPrCounter}.ccgSmooth, 'r');
                    line([0 0], ylim, 'Color', 'r');
                    line([ccgPars{validCntrPrCounter}.firstPeak, ccgPars{validCntrPrCounter}.firstPeak], ylim, 'Color', 'c');
                    %subplot(224);
                    reportfig(gcf, [mfilename], 0, ['filebase : ', gt.filebase, '  trial name : ', gt.trialName, '   ', roi]);
                    clf
                    %%%%%%
                    % offset as afunction of rate pk distance
                    dvsoff = [dvsoff; p(validCntrPrCounter), ccgPars{validCntrPrCounter}.firstPeak];
                end
            end
            out.ccg = ccgPars;
            out.dvsoff = dvsoff;
        end
        save(['~/data/analysis/', datasetType, '/', filebase, '/' filebase, '.', mfilename, '.mat'], 'out');
        
      case 'load'
        load(['~/data/analysis/', datasetType, '/', filebase, '/' filebase, '.', mfilename, '.mat'], 'out');
    end
end

