function out = MultiPeakPFDistance(gt, varargin)
% out = MultiPeakPFDistance(gpf, occupancy, varargin)
% 
    if isempty(gt.pfObject), gt = gt.LoadPF; end
    gpf = gt.pfObject;
    defClu = gpf.acceptedUnits;
    defCellPairs = nchoosek(defClu, 2);
    [cellPairs, IF_SMRM, IF_COMPUTE_CCG, nSTD, areaThreshFactor, occThreshFac, state, binSize, maxTimeLag] = ...
        DefaultArgs(varargin, {defCellPairs, 1, 0, 3, 0.5, 2, 'RUN', 10e-3, 200e-3});
    
    cellIds = unique(cellPairs(:));
    nCells = length(cellIds);
    smthRateMaps = gpf.smoothRateMap(:, :, nCells);
    occupancy = Occupancy(gt);
    occThreshold = occThreshFac * std(occupancy(:));
    occupancy(occupancy <= occThreshold) = 0;
    validXYBins = Ind2Sub(size(occupancy), find(occupancy > 0));
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
                SELECTED_CELL(kCell) = false;
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
        cellIds = [];
    else
        cellIds = cellIds(SELECTED_CELL);
        out.cntrVertices = cntrVertices;
        out.cntrPeaks = cntrPeaks;
    end
    %%%%%%  CCG 
    if IF_COMPUTE_CCG
        if length(cellIds) > 1, cellPairs = nchoosek(cellIds, 2); 
        else, fprintf('\n no cells with adequate sampling found'); return; end
        [res, residx] =SelectPeriods(gt.res, ConvertFs(gt.trialPeriods, gt.lfpSampleRate, gt.sampleRate), 'd', 1, 1);
        clu = gt.clu(residx);
        res = ConvertFs(res, gt.sampleRate, gt.trackingSampleRate);
        clu = clu(ismember(clu, cellIds));
        res = res(ismember(clu, cellIds));
        binnedPos = BinPos(gt);
        halfBins = round(maxTimeLag / binSize); % number of bins on each side of zero lag
        binSize = round(binSize * gt.sampleRate);
        nCellPairs = size(cellPairs, 1);
        stsp = gt.LoadStatePeriods([],gt.trackingSampleRate, 1);
        stsp = stsp - ConvertFs(gt.trialPeriods(1), gt.lfpSampleRate, gt.trackingSampleRate);
        counterPk = 0;
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
                
                for kCntrPr = 1 : size(cntrPairs, 1)
                    kCntrPr 
                    % OVERLAP = any(InPolygon(cntrA{cntrPairs(kCntrPr, 1)}, cntrB{cntrPairs(kCntrPr, 2)}));
                    % if the sub cntrs overlap
                    pkDistAB(counterPk + 1) = norm( pkA(cntrPairs(kCntrPr, 1), :) - pkB(cntrPairs(kCntrPr, 1), :));
                    if any(InPolygon(cntrA{cntrPairs(kCntrPr, 1)}, cntrB{cntrPairs(kCntrPr, 2)})); 
                        POS_IN_CNTR_A = InPolygon(binnedPos, cntrA{cntrPairs(kCntrPr, 1)});
                        POS_IN_CNTR_B = InPolygon(binnedPos, cntrB{cntrPairs(kCntrPr, 2)});
                        POS_IN_AB = POS_IN_CNTR_A | POS_IN_CNTR_B;
                        inout = InOut(POS_IN_AB);
                        timesInAB = inout{1};
                        timesInAB = IntersectRanges(timesInAB, stsp);
                        res1 = SelectPeriods(res(clu == cellPairs(mCellPair, 1)), timesInAB, 'd');
                        res2 = SelectPeriods(res(clu == cellPairs(mCellPair, 2)), timesInAB, 'd');
                        if length(res1) > 5 & length(res2) > 5
                            mClu = [ones(length(res1), 1) * cellPairs(mCellPair, 1); ones(length(res2), 1) * cellPairs(mCellPair, 2)];
                            [mCCg, ccgTAx, ~] = myCCG([res1; res2], mClu, binSize, halfBins, gt.sampleRate, cellPairs(mCellPair, :), 'count', [], 1);
                            xsc%%%% DISPLAY %%%%%
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
                            plot(pkB(cntrPairs(kCntrPr, 2), 1), pkB(cntrPairs(kCntrPr, 2), 2), 'g*')
                            subplot(2, 2, 1); hold on;
                            line([0 0], ylim, 'Color', 'r');
                            %%%%%%
                            keyboard;
                            clf;                   
                        end
                    end
                end
            else
                PF_OVERLAP(mCellPair) = false;
            end
        end
        keyboard;
    end
end
