function out = CatTrialCCG(filebase, roi,  arena, state, varargin)

     [cellPairs, datasetType] = DefaultArgs(varargin, {[], 'kenji'});

    [matches, fbRoi] =  SearchKenji(filebase);
    if ~any(cellfun(@strcmp, fbRoi, repmat({roi}, size(fbRoi)))); return; end
    trialNames = matches(cellfun(@strcmp, matches(:, 3), repmat({arena}, size(matches(:, 3)))), 2);
    trialPeriods = [];
    load(['~/data/analysis/', datasetType, '/Contours', GenFiletag(roi, arena), 'mat'], 'cntrs', 'filebases');
    [stableCntrs, cellIds] = StableCntrs(cntrs(cellfun(@strcmp, filebases, repmat({filebase}, size(filebases))), :), roi, arena);
    cellIds = cell2mat(cellIds);
    for mTr = 1 : length(trialNames)
        gt = GenericTrial(filebase, trialNames{mTr});
        if length(cellIds) > 1, cellPairs = nchoosek(cellIds, 2); else, fprintf('\n no pairs found'); return; end
        mpf = stableCntrs{mTr};
        cntrVertices = mpf.cntrVertices;
        cntrPeaks = mpf.cntrPeaks;
        trialPeriods = [trialPeriods; gt.LoadStatePeriods(state, 0, 0)];
        [res, clu] = gt.LoadStateRes(state, 1);
        clu = clu(ismember(clu, cellIds));
        res = res(ismember(clu, cellIds));
        binnedPos = BinPos(gt);
        nCellPairs = size(cellPairs, 1);
        stsp = gt.LoadStatePeriods(state, gt.lfpSampleRate, 1);
        stsp = stsp - gt.trialPeriods(1);
        counterPk = 0;
        PF_OVERLAP = false(nCellPairs, 1);
        ccgPars = [];
        ncntrpr(mTr) = 0;
        catRes{1} = [];
        catClu{1} = [];
        cntrprCounter = 0;
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
                ncntrpr(mTr) = ncntrpr(mTr) + size(cntrPairs, 1);
                tcatRes = [];
                tcatClu = [];

                for kCntrPr = 1 : size(cntrPairs, 1)
                    kCntrPr 
                    % if the sub cntrs overlap
                    %pkDistAB(validCntrCnt + 1) = norm( pkA(cntrPairs(kCntrPr, 1), :) - pkB(cntrPairs(kCntrPr, 2), :));
                    if any(InPolygon(cntrA{cntrPairs(kCntrPr, 1)}, cntrB{cntrPairs(kCntrPr, 2)})); 
                        POS_IN_CNTR_A = InPolygon(binnedPos, cntrA{cntrPairs(kCntrPr, 1)});
                        POS_IN_CNTR_B = InPolygon(binnedPos, cntrB{cntrPairs(kCntrPr, 2)});
                        POS_IN_AB = POS_IN_CNTR_A | POS_IN_CNTR_B;
                        inout = InOut(POS_IN_AB);
                        timesInAB = inout{1};
                        timesInAB = ConvertFs(IntersectRanges(timesInAB, ConvertFs(stsp, gt.lfpSampleRate, gt.trackingSampleRate)), gt.trackingSampleRate, gt.sampleRate);
                        res1 = SelectPeriods(res(clu == cellPairs(mCellPair, 1)), timesInAB, 'd', 1, 1);
                        res2 = SelectPeriods(res(clu == cellPairs(mCellPair, 2)), timesInAB, 'd', 1, 1);
                        mClu = [ones(length(res1), 1) * cellPairs(mCellPair, 1); ones(length(res2), 1) * cellPairs(mCellPair, 2)];
                        cntrprCounter = cntrprCounter + 1;
                        if mTr == 1, catRes{cntrprCounter} = []; catClu{cntrprCounter} = []; end
                        catRes{cntrprCounter} = [catRes{cntrprCounter}; res1; res2];
                        catClu{cntrprCounter} = [catClu{cntrprCounter}; mClu];
                        pkDistAB(cntrprCounter) = norm( pkA(cntrPairs(kCntrPr, 1), :) - pkB(cntrPairs(kCntrPr, 2), :));
                        PF_OVERLAP(mCellPair) = true;
                    end
                end
            end
        end
    end
keyboard;
end

        

%                         if length(res1) > 10 & length(res2) > 10
%                             PF_OVERLAP(mCellPair) = true;
%                             mClu = [ones(length(res1), 1) * cellPairs(mCellPair, 1); ones(length(res2), 1) * cellPairs(mCellPair, 2)];
%                             options = struct('type', 'jitter', 'winSize', 20e-3, 'nResamples', 1e2);
%                             h1 = figure(79);
                            
%                             tccg = CCGPars([res1; res2], mClu, gt.sampleRate, [], options, binSize, maxTimeLag, [],  0, IF_PLOT);
%                             validCntrCnt = validCntrCnt + 1;
%                             ccgPars{mCellPair, validCntrCnt} = tccg{1};
%                             %%%% DISPLAY %%%%%
%                             if IF_PLOT
%                                 cellPairs(mCellPair, :)
%                                 figure(h1);
%                                 subplot(2,2,4);
%                                 cla;
%                                 imagesc(occupancy); colormap('gray');
%                                 set(gca, 'ydir', 'normal');
%                                 hold on;
%                                 plot(cntrA{cntrPairs(kCntrPr, 1)}(:,1), cntrA{cntrPairs(kCntrPr, 1)}(:,2),'r.-');
%                                 hold on
%                                 plot(pkA(cntrPairs(kCntrPr, 1), 1), pkA(cntrPairs(kCntrPr, 1), 2), 'g*')
%                                 xlim([1 50]);
%                                 ylim([1 50]);
%                                 plot(cntrB{cntrPairs(kCntrPr, 2)}(:,1), cntrB{cntrPairs(kCntrPr, 2)}(:,2),'y.-');
%                                 hold on
%                                 plot(pkB(cntrPairs(kCntrPr, 2), 1), pkB(cntrPairs(kCntrPr, 2), 2), 'g*');
%                                 subplot(2, 2, 1); hold on;
%                                 plot(ccgPars{mCellPair, validCntrCnt}.smoothTAx, ccgPars{mCellPair, validCntrCnt}.ccgSmooth, 'r');
%                                 line([0 0], ylim, 'Color', 'r');
%                                 line([ccgPars{mCellPair, validCntrCnt}.firstPeak, ccgPars{mCellPair, validCntrCnt}.firstPeak], ylim, 'Color', 'c');
%                                 reportfig(h1, [mfilename], 0, [gt.filebase, '--' gt.trialName]);
%                                 clf
%                                 %%%%%%
%                                 figure(37)
%                                 % offset as afunction of rate pk distance
%                                 dvsoff = [dvsoff; pkDistAB(validCntrCnt), ccgPars{mCellPair, validCntrCnt}.firstPeak];
%                                 plot(pkDistAB(validCntrCnt), ccgPars{mCellPair, validCntrCnt}.firstPeak, '*');
%                                 hold on;
%                             end
%                         end
%                     end
%                 end
                
%             else
%                disp('here');
%             end
%             %}
    


   