function [stableCntrs, selectedClu, allFbs, processedFbs] = StableCntrs(inCntrs, roi, arena)
% StableCntrs(contours, roi, arena)
% contours is the cell array returned by GetContours
% this scripts looks for stable sub-fields across trials in similar arena, if only one trial exists, 
% returns cntrs identical to  inCntrs

%load('~/thesis/temp/poolOffset_cntrs');
    sk.roi = roi;
    sk.arena = arena;
    matches = SearchKenji(sk);
    filebases = unique(matches(:, 1));

    passedCntrs = inCntrs;
    inCntrs = inCntrs(logical(sum(~cellfun(@isempty, inCntrs), 2)), :); % remove empty rows
    selectedFBs = filebases(logical(sum(~cellfun(@isempty, inCntrs), 2)));
    nFbs = length(selectedFBs);
    allFbs = selectedFBs;
    nTrs = sum(~cellfun(@isempty, inCntrs), 2);
    inCntrs = inCntrs(nTrs > 1, :);
    selectedFBs = selectedFBs(nTrs > 1);
    processedFbs = selectedFBs;
    identicalFbIdx = find(nTrs ==1); % fbs with only 1 trial
    diffFbIdx = find(nTrs > 1); % fbs with nore than 1 trial
    maxTrs = max(nTrs);
    nTrs = nTrs(nTrs > 1);
    % find pks of cntrs which are stable acorss sessions
    for mBase = 1 : size(inCntrs, 1)
        DetectCmnPks = GenDetectCmnPks(2, 'ismemberf', 3);
        cntrPeakNames = genvarname(cellstr(repmat('cntrPeak', nTrs(mBase), 1)));
        load(['~/data/analysis/kenji/', filebases{mBase}, '/', filebases{mBase}, GenFiletag(roi, arena), 'commonClus.mat']);
        trPairs = nchoosek(1 : nTrs(mBase), 2);
        trPairs(trPairs(:, 1) ~= 1 ,:) = [] ;
        for lTr = 1 : nTrs(mBase)
            eval([cntrPeakNames{lTr} '=inCntrs{mBase, lTr}.cntrPeaks;']);
        end
        for kTrPr = 1 : size(trPairs, 1)
            evalStr = [];
            for kPr = 1 : 2
                evalStr = [evalStr, ',' cntrPeakNames{trPairs(kTrPr, kPr)}];
                IS_EMPTY = eval(['isempty(' cntrPeakNames{trPairs(kTrPr, kPr)} ');']);
            end
            if IS_EMPTY
                kTF{kTrPr} = {[]}; kLoc{kTrPr} = {[]};
            else
                [kTF{kTrPr}, kLoc{kTrPr}] = eval(['cellfun(DetectCmnPks' evalStr, ',', '''uniformoutput''' ', 0);']);
            end
        end
        nTrPairs = size(trPairs, 1);
        if nTrPairs > 1
            for mTrPr = 1 : nTrPairs - 1
                if mTrPr == 1
                    mTF{mBase} = cellfun(@and, kTF{mTrPr}, kTF{mTrPr + 1}, 'UNIFORMOUTPUT', 0); 
                else
                    mTF{mBase} = cellfun(@and, mTF{mBase}, kTF{mTrPr + 1},'UNIFORMOUTPUT', 0);
                end
            end
        else
            mTF{mBase} = kTF{1};
        end
        STABLE_CELLS{mBase}  = logical(cellfun(@sum, mTF{mBase}));
        nStableCells = sum(STABLE_CELLS{mBase});
        nStableCntrs = sum(cellfun(@sum, mTF{mBase}));
        for mTrPr = 1 : nTrPairs
            for lPr = 1 : 2
                idx = [];
                if lPr== 1
                    idx = mTF{mBase}(STABLE_CELLS{mBase});
                else
                    tempKloc = kLoc{mTrPr}(STABLE_CELLS{mBase});
                    tempLkIdx = mTF{mBase}(STABLE_CELLS{mBase});
                    for lCell = 1 : nStableCells
                        idx{lCell} = tempKloc{lCell}(tempLkIdx{lCell});
                    end
                end
                cntrVertices = inCntrs{mBase, trPairs(mTrPr, lPr)}.cntrVertices(STABLE_CELLS{mBase});
                cntrPeaks = inCntrs{mBase, trPairs(mTrPr, lPr)}.cntrPeaks(STABLE_CELLS{mBase});
                cntrs = struct('cntrVertices', [], 'cntrPeaks', []);
                for kCell = 1 : nStableCells
                    kIdx = idx{kCell};
                    if ~isempty(kIdx)
                        if lPr ~= 1, kIdx(kIdx == 0) = []; end
                        tempCntr = cntrVertices{kCell};
                        tempPeaks = cntrPeaks{kCell};
                        cntrs.cntrVertices{kCell} = tempCntr(kIdx);
                        cntrs.cntrPeaks{kCell} = tempPeaks(kIdx, :);
                        %                    else
                        %cntrs = struct('cntrVertices', [], 'cntrPeaks', []);
                    end
                end
            cmnCntrs{mBase, trPairs(mTrPr, lPr)} = cntrs;
        end
    end
    clu{mBase} = commonClus(STABLE_CELLS{mBase});
    end
    stableCntrs = cell(nFbs, maxTrs);
    stableCntrs(identicalFbIdx) = passedCntrs(identicalFbIdx);
    stableCntrs(diffFbIdx, :) = cmnCntrs; 
    selectedClu = cell(1, nFbs);
    selectedClu(diffFbIdx) = clu;

    
     
    
    
    
    
% %% Rate remapping
% for kBase = 1 : length(selectedFBs)
%     kTrials = matches(strcmp(matches(:, 1), selectedFBs{kBase}), 2);
%     kClu = y.clu{kBase};
%     refTr = 1;
%     for lTr = 1 : length(kTrials)
%         gt = GenericTrial(selectedFBs{kBase}, kTrials{lTr});
%         gt = gt.LoadPF;
%         srm = reshape(gt.pfObject.smoothRateMap(:,:, ismember(gt.pfObject.acceptedUnits, y.clu{kBase})), [], length(kClu));
%         for mCell = 1 : length(kClu);
%             mCellIdx = Sub2Ind([length(gt.pfObject.xBin), length(gt.pfObject.yBin)], y.cmnCntrs{kBase, lTr}.cntrPeaks{mCell});
%             kSrm = srm(:, mCell);
%             mCellPk{mCell} = kSrm(mCellIdx);
%             tempPk(:, lTr) = mCellPk{mCell};
%         end
%         %        rmpk{lTr} = mCellPk;

% %         if lTr == refTr
% %             refTrRMPk = rmpk{lTr};
% %         end
%     end
%     for mTr = 1: length(kTrials)
%         if mTr ~= refTr
%             plot(tempPk(:, refTr), tempPk(:, mTr), '.');
%             hold all;
%            end

% end

%% compute distance between peaks of stable cntrs

% for kBase = 1 : size(inCntrs, 1)
%     for kTr = 1 : nTrs(kBase)
%         for mCells  = 1 : length(cmnCntrs{kBase, kTr}.cntrVertices)
%             sc1 = c1{ii};
%             sc2 = c2{ii};sc3 = c3{ii};
%             for jj = 1 : length(sc1)
%                 plot(sc1{jj}(:,1), sc1{jj}(:,2));
%                 hold on;
%                 plot(sc2{jj}(:,1), sc2{jj}(:,2), 'r');plot(sc3{jj}(:,1), sc3{jj}(:,2), 'g');
%             end
%             waitforbuttonpress
%         end
%     end
% end

end