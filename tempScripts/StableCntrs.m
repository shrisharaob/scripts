
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






    out = out(logical(sum(~cellfun(@isempty, out), 2)), :); % remove empty rows
                                                            %    selectedFBs = filebases(logical(sum(~cellfun(@isempty, out), 2)));
    nTrs = sum(~cellfun(@isempty, out), 2);
    out = out(nTrs > 1, :);
    nTrs = nTrs(nTrs > 1);
    % find pks of cntrs which are stable acorss sessions
    for mBase = 1 : size(out, 1)
        DetectCmnPks = GenDetectCmnPks(2, 'ismemberf', 3);
        cntrPeakNames = genvarname(cellstr(repmat('cntrPeak', nTrs(mBase), 1)));

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
                        keyboard;
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
    end

    