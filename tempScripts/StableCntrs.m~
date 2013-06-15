
%     for kBase = 1 : length(filebases)
%         trialNames = TrialNames(filebases{kBase}, datasetType, roi, arena);
%         fprintf([filebases{kBase} '\n']);
%         for mTr = 1 : length(trialNames)
%             fprintf(['\t' trialNames{mTr} '\n']);
%             gt =  GenericTrial(filebases{kBase}, trialNames{mTr});
%             try
%                 % load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(arena, roi), 'CCG.mat' ]);
% %                 load(['~/data/analysis/kenji/', filebases{kBase}, '/', filebases{kBase}, '.', trialNames{mTr}, GenFiletag(arena, roi), 'CCG.mat' ]);
% %                 offset{kBase, mTr} = abs(fout.offset);
%                 commonClus = gt.LoadCommonClus;
%                 out{kBase, mTr} = gt.MultiPeakPFDistance(nchoosek(commonClus, 2));
%                 trNo = trNo + 1;
%                 trNames{trNo} = trialNames{mTr}; 
%             catch err
%                 fprintf(['\n error ' trialNames{mTr}, '\n']);
%             end     
%             %kBase 
%         end
%     end
load('~/thesis/temp/poolOffset_cntrs');
matches = SearchKenji;
filebases = unique(matches(:, 1));
out = out(logical(sum(~cellfun(@isempty, out), 2)), :); % remove empty rows
selectedFBs = filebases(logical(sum(~cellfun(@isempty, out), 2)));
nTrs = sum(~cellfun(@isempty, out), 2);
out = out(nTrs > 1, :);
selectedFBs = selectedFBs(nTrs > 1);
nTrs = nTrs(nTrs > 1);
% find pks of cntrs which are stable acorss sessions
for mBase = 1 : size(out, 1)
    DetectCmnPks = GenDetectCmnPks(2, 'ismemberf', 3);
    cntrPeakNames = genvarname(cellstr(repmat('cntrPeak', nTrs(mBase), 1)));
    load(['~/data/analysis/kenji/', filebases{mBase}, '/', filebases{mBase}, '.CA3.bigSquare.commonClus.mat']);
    trPairs = nchoosek(1 : nTrs(mBase), 2);
    trPairs(trPairs(:, 1) ~= 1 ,:) = [] ;
    for lTr = 1 : nTrs(mBase)
        eval([cntrPeakNames{lTr} '=out{mBase, lTr}.cntrPeals;']);
    end
    for kTrPr = 1 : size(trPairs, 1)
        evalStr = [];
        for kPr = 1 : 2
            evalStr = [evalStr, ',' cntrPeakNames{trPairs(kTrPr, kPr)}];
        end
        [kTF{kTrPr}, kLoc{kTrPr}] = eval(['cellfun(DetectCmnPks' evalStr, ',', '''uniformoutput''' ', 0);']);
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
            cntrVertices = out{mBase, trPairs(mTrPr, lPr)}.cntrVertices(STABLE_CELLS{mBase});
            cntrPeaks = out{mBase, trPairs(mTrPr, lPr)}.cntrPeals(STABLE_CELLS{mBase});
            for kCell = 1 : nStableCells
                kIdx = idx{kCell};
                if ~isempty(kIdx)
                    if lPr ~= 1, kIdx(kIdx == 0) = []; end
                    tempCntr = cntrVertices{kCell};
                    tempPeaks = cntrPeaks{kCell};
                    cntrs.cntrVertices{kCell} = tempCntr(kIdx);
                    cntrs.cntrPeaks{kCell} = tempPeaks(kIdx, :);
                else
                    cntrs = {};
                end
            end
            cmnCntrs{mBase, trPairs(mTrPr, lPr)} = cntrs;
        end
    end
clu{mBase} = commonClus(STABLE_CELLS{mBase});
end
y.cmnCntrs = cmnCntrs;
y.clu = clu;

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

% for kBase = 1 : size(out, 1)
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

