function out = MultiPeakPFDistance(gt, roi, arena, varargin)
% out = MultiPeakPFDistance(gpf, occupancy, varargin)
% [cellPairs, IF_SMRM, IF_COMPUTE_CCG, nSTD, areaThreshFactor, occThreshFac, state, binSize, maxTimeLag, IF_PLOT] = ...

    if isempty(gt.pfObject), gt.LoadPF; end
    gpf = gt.pfObject;
    %    defClu = gpf.acceptedUnits;
    defClu = gt.LoadCommonClus(roi, arena);
    out.pkDist = [];
    if isempty(defClu), return; end
    if length(defClu) < 2, return; end
    defCellPairs = nchoosek(defClu, 2);
    [cellPairs, IF_SMRM, IF_COMPUTE_CCG, IF_CAT_TRIALS, nSTD, areaThreshFactor, occThreshFac, state, binSize, maxTimeLag, IF_PLOT] = ...
        DefaultArgs(varargin, {defCellPairs, 1, 0, 0, 3, 0.5, 0, 'RUN', 20e-3, 1000e-3, true});
    
    cellIds = unique(cellPairs(:));
    nCells = length(cellIds);
    smthRateMaps = gpf.smoothRateMap(:, :, nCells);
    occupancy = Occupancy(gt);
    occThreshold = occThreshFac * std(occupancy(:));
    occupancy(occupancy <= occThreshold) = 0;
    validXYBins = Ind2Sub(size(occupancy), find(occupancy > 0));
    SELECTED_CELL = false(1, nCells);
    dvsoff = [];
    % detect individual contours and their centers satisfying minimum area criterion
    for kCell = 1 : nCells
        if IF_SMRM
            idx = ismember(gpf.acceptedUnits, cellIds(kCell));
            if ~isempty(idx) & sum(idx) ~= 0
                kSmoothRM = gpf.smoothRateMap(:, :, idx);
                rateThresh = nSTD * std(kSmoothRM(:));
                kPk = LocalMinima2(-1 * kSmoothRM', -1 * rateThresh, 2);
                kContr = contour(kSmoothRM, [1, 1] .* rateThresh);
                [nRows, nClmns] = size(kContr);
                nContrs = 0;
                IS_DONE = 0;
                nVals = 0;
                while ~IS_DONE
                    nContrs = nContrs + 1;
                    mClm(nContrs) = 1 + sum(nVals) + (nContrs - 1);
                    if mClm(nContrs) > nClmns, 
                        mClm(nContrs) = nClmns;
                        nContrs = nContrs - 1;
                        IS_DONE = 1; continue; 
                    end
                    nVals(nContrs) = kContr(2, mClm(nContrs));
                end
                for mCntr = 1 : nContrs
                    tempCntr = kContr(:, mClm(mCntr) + 1 : mClm(mCntr + 1) - 1 * ~(mCntr == nContrs));
                    % if the countours are not closed, complete it
                    if ~all(tempCntr(:, 1) == tempCntr(:, end)), tempCntr(:, end + 1) = tempCntr(:, 1); end 
                    mCntrVertices{mCntr} = tempCntr';
                    mArea(mCntr) = polyarea(mCntrVertices{mCntr}(:, 1), mCntrVertices{mCntr}(:, 2));
                    mPkInCntr = find(InPolygon(kPk, mCntrVertices{mCntr})); % pks inside the cntrs
                    if length(mPkInCntr) > 1
                        pkIdx = Sub2Ind(size(kSmoothRM), kPk(mPkInCntr, :));
                        [~, maxPkIdx] = max(kSmoothRM(pkIdx));
                        mPkInCntr = mPkInCntr(maxPkIdx, :);  
                    end
                    if isempty(mPkInCntr)
                        cntrPk(mCntr, :) = [nan, nan];
                    else
                        cntrPk(mCntr, :) = kPk(mPkInCntr, :);
                    end
                end
                [maxArea, maxCntrId] = max(mArea);
                areaThresh = areaThreshFactor * maxArea; % discard pf subfield if area less than thresh
                IS_VALID_CNTR{kCell} = mArea >= areaThresh & ~isnan(cntrPk(1,1));
                cntrVertices{kCell} = mCntrVertices(IS_VALID_CNTR{kCell});
                cntrPeaks{kCell} = cntrPk(IS_VALID_CNTR{kCell}, :);
                for lValidCntr =  1 : sum(IS_VALID_CNTR{kCell})
                    % the unit is selected if it has more than thresh occupancy
                    SELECTED_CELL(kCell) = SELECTED_CELL(kCell) | any(InPolygon(validXYBins, mCntrVertices{1}));
                end
                clear mArea mCntrVertices cntrPk
                threshMask(:, :, kCell) = kSmoothRM > rateThresh;
            end
        end
    end
    if ~(exist('SELECTED_CELL', 'var')), 
        out.cntrVertices = [];
        out.cntrPeaks = [];
        out.cluId = [];
        cellIds = [];
    else
        out.cntrVertices = cntrVertices(SELECTED_CELL);
        out.cntrPeaks = cntrPeaks(SELECTED_CELL);
        out.cluId = cellIds(SELECTED_CELL);
    end
   
    %% Peak distances
    nCellPairs = size(cellPairs, 1);
    for mCellPair = 1 : nCellPairs
        validCntrCnt = 0;
        cntrA = cntrVertices{cellIds == cellPairs(mCellPair, 1)}; 
        cntrB = cntrVertices{cellIds == cellPairs(mCellPair, 2)};
        nCntrA = length(cntrA);
        nCntrB = length(cntrB);
        mCellPair
        if nCntrA >= 1 && nCntrB >= 1, 
            cntrPairs = nchoosek([1 : nCntrA, 1 : nCntrB], 2); % all pairs of selected sub contours
            cntrPairs(cntrPairs(:, 1) > nCntrA, :) = [];
            cntrPairs(cntrPairs(:, 2) > nCntrB, :) = [];
            cntrPairs = sortrows(unique(cntrPairs, 'rows')); 
            pkA = cntrPeaks{cellIds == cellPairs(mCellPair, 1)};
            pkB = cntrPeaks{cellIds == cellPairs(mCellPair, 2)};
            for kCntrPr = 1 : size(cntrPairs, 1)
                kCntrPr 
                validCntrCnt = validCntrCnt + 1;
                pkDistAB(validCntrCnt) = norm( pkA(cntrPairs(kCntrPr, 1), :) - pkB(cntrPairs(kCntrPr, 2), :));
                selectedCellpairs(validCntrCnt, :) = cellPairs(mCellPair, :);
                pkAB(validCntrCnt, :) = [pkA(cntrPairs(kCntrPr, 1), :) , pkB(cntrPairs(kCntrPr, 2), :)];
            end
        else
            selectedCellpairs = []; pkDistAB = []; pkAB = [];
        end
    end

    out.pkDist = [selectedCellpairs, pkDistAB', pkAB];
    %%  CCG for spikes fired when the animal is in the region covered by contour pairs
    if IF_COMPUTE_CCG
        if length(cellIds) > 1, cellPairs = nchoosek(cellIds, 2); 
        else, fprintf('\n no cells with adequate sampling found'); return; end
        if IF_CAT_TRIALS
            if isempty(gt.clu), gt.LoadCR; end
            catTrPeriods = CatTrials(gt.filebase, roi, arena, state, gt.sampleRate);
            [res, resIdx] = SelectPeriods(gt.res, catTrPeriods, 'd', 1, 1);
            clu  = gt.clu(resIdx);
        else
            [res, clu] = gt.LoadStateRes(state, 1);
        end
        clu = clu(ismember(clu, cellIds));
        res = res(ismember(clu, cellIds));
        binnedPos = BinPos(gt);
        nCellPairs = size(cellPairs, 1);
        stsp = gt.LoadStatePeriods(state, gt.lfpSampleRate, 1);
        stsp = stsp - gt.trialPeriods(1);
        counterPk = 0;
        PF_OVERLAP = false(nCellPairs, 1);
        ccgPars = [];
        for mCellPair = 1 : nCellPairs
            cntrA = cntrVertices{cellIds == cellPairs(mCellPair, 1)}; 
            cntrB = cntrVertices{cellIds == cellPairs(mCellPair, 2)};
            nCntrA = length(cntrA);
            nCntrB = length(cntrB);
            mCellPair
            if nCntrA >= 1 && nCntrB >= 1, 
                cntrPairs = nchoosek([1 : nCntrA, 1 : nCntrB], 2); % all pairs of selected sub contours
                cntrPairs(cntrPairs(:, 1) > nCntrA, :) = [];
                cntrPairs(cntrPairs(:, 2) > nCntrB, :) = [];
                cntrPairs = sortrows(unique(cntrPairs, 'rows')); 
                pkA = cntrPeaks{cellIds == cellPairs(mCellPair, 1)};
                pkB = cntrPeaks{cellIds == cellPairs(mCellPair, 2)};
                validCntrCnt = 0;
                for kCntrPr = 1 : size(cntrPairs, 1)
                    kCntrPr 
                    % if the sub cntrs overlap
                    pkDistAB(validCntrCnt + 1) = norm( pkA(cntrPairs(kCntrPr, 1), :) - pkB(cntrPairs(kCntrPr, 2), :));
                    if any(InPolygon(cntrA{cntrPairs(kCntrPr, 1)}, cntrB{cntrPairs(kCntrPr, 2)})); 
                        POS_IN_CNTR_A = InPolygon(binnedPos, cntrA{cntrPairs(kCntrPr, 1)});
                        POS_IN_CNTR_B = InPolygon(binnedPos, cntrB{cntrPairs(kCntrPr, 2)});
                        POS_IN_AB = POS_IN_CNTR_A | POS_IN_CNTR_B;
                        inout = InOut(POS_IN_AB);
                        timesInAB = inout{1};
                        timesInAB = ConvertFs(IntersectRanges(timesInAB, ConvertFs(stsp, gt.lfpSampleRate, gt.trackingSampleRate)), gt.trackingSampleRate, gt.sampleRate);
                        res1 = SelectPeriods(res(clu == cellPairs(mCellPair, 1)), timesInAB, 'd');
                        res2 = SelectPeriods(res(clu == cellPairs(mCellPair, 2)), timesInAB, 'd');
                        if length(res1) > 10 & length(res2) > 10
                            PF_OVERLAP(mCellPair) = true;
                            mClu = [ones(length(res1), 1) * cellPairs(mCellPair, 1); ones(length(res2), 1) * cellPairs(mCellPair, 2)];
                            options = struct('type', 'jitter', 'winSize', 20e-3, 'nResamples', 1e2);
                            h1 = figure(79);
                            tccg = CCGPars([res1; res2], mClu, gt.sampleRate, [], options, binSize, maxTimeLag, [],  0, IF_PLOT);
                            validCntrCnt = validCntrCnt + 1;
                            ccgPars{mCellPair, validCntrCnt} = tccg{1};
                            %%%% DISPLAY %%%%%
                            if IF_PLOT
                                cellPairs(mCellPair, :)
                                figure(h1);
                                subplot(2,2,4);
                                cla;
                                imagesc(occupancy); colormap('gray');
                                set(gca, 'ydir', 'normal');
                                hold on;
                                plot(cntrA{cntrPairs(kCntrPr, 1)}(:,1), cntrA{cntrPairs(kCntrPr, 1)}(:,2),'r.-');
                                hold on
                                plot(pkA(cntrPairs(kCntrPr, 1), 1), pkA(cntrPairs(kCntrPr, 1), 2), 'g*')
                                xlim([1 50]);
                                ylim([1 50]);
                                plot(cntrB{cntrPairs(kCntrPr, 2)}(:,1), cntrB{cntrPairs(kCntrPr, 2)}(:,2),'y.-');
                                hold on
                                plot(pkB(cntrPairs(kCntrPr, 2), 1), pkB(cntrPairs(kCntrPr, 2), 2), 'g*');
                                subplot(2, 2, 1); hold on;
                                plot(ccgPars{mCellPair, validCntrCnt}.smoothTAx, ccgPars{mCellPair, validCntrCnt}.ccgSmooth, 'r');
                                line([0 0], ylim, 'Color', 'r');
                                line([ccgPars{mCellPair, validCntrCnt}.firstPeak, ccgPars{mCellPair, validCntrCnt}.firstPeak], ylim, 'Color', 'c');
                                reportfig(h1, [mfilename], 0, [gt.filebase, '--' gt.trialName]);
                                clf
                                %%%%%%
                                figure(37)
                                % offset as afunction of rate pk distance
                                dvsoff = [dvsoff; pkDistAB(validCntrCnt), ccgPars{mCellPair, validCntrCnt}.firstPeak];
                                plot(pkDistAB(validCntrCnt), ccgPars{mCellPair, validCntrCnt}.firstPeak, '*');
                                hold on;
                            end
                        end
                    end
                end
                
            else
               disp('here');
            end
        end
        out.ccgPars = ccgPars;
        out.PF_OVERLAP = PF_OVERLAP;
        out.dvsoff = dvsoff;
    end
end
