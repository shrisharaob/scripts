function PoolOffset(varargin)

    [datasetType, roi, arena, nResample] = DefaultArgs(varargin, {'kenji', 'CA3', 'bigSquare', 1e2});
    switch datasetType
        case 'kenji'
          sK.roi = roi;
          sK.arena = arena;
          filebases = SearchKenji(sK);
          filebases = unique(filebases(:, 1));
    end
    trNo = 0;
    for kBase = 1 : length(filebases)
        trialNames = TrialNames(filebases{kBase}, datasetType, roi, arena);
        fprintf([filebases{kBase} '\n']);
        for mTr = 1 : length(trialNames)
            fprintf(['\t' trialNames{mTr} '\n']);
            gt =  GenericTrial(filebases{kBase}, trialNames{mTr});
            try
                % load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(arena, roi), 'CCG.mat' ]);
%                 load(['~/data/analysis/kenji/', filebases{kBase}, '/', filebases{kBase}, '.', trialNames{mTr}, GenFiletag(arena, roi), 'CCG.mat' ]);
%                 offset{kBase, mTr} = abs(fout.offset);
                commonClus = gt.LoadCommonClus(roi, arena);
                if  isempty(commonClus)
                    cntrs{kBase, mTr} = struct('cntrVertices', [], 'cntrPeaks', []);
                elseif length(commonClus) < 2
                    cntrs{kBase, mTr} = struct('cntrVertices', [], 'cntrPeaks', []);
                else
                    cntrs{kBase, mTr} = gt.MultiPeakPFDistance(nchoosek(commonClus, 2));
                end
                trNo = trNo + 1;
                trNames{trNo} = trialNames{mTr}; 
            catch err
                fprintf(['\n error ' trialNames{mTr}, '\n']);
                keyboard;
            end     
        end
    end
    save(['~/data/analysis/', datasetType, '/Contours', GenFiletag(roi, arena), 'mat'], 'cntrs');
    keyboard;
%     out = out(logical(sum(~cellfun(@isempty, out), 2)), :); % remove empty rows
%     selectedFBs = filebases(logical(sum(~cellfun(@isempty, out), 2)));
%     nTrs = sum(~cellfun(@isempty, out), 2);
%     out = out(nTrs > 1, :);
%     nTrs = nTrs(nTrs > 1);
%     % find pks of cntrs which are stable acorss sessions
%     for mBase = 1 : size(out, 1)
%         DetectCmnPks = GenDetectCmnPks(nTrs(mBase), 'ismemberf', 3);
%         cntrPeakNames = genvarname(cellstr(repmat('cntrPeak', nTrs(mBase), 1)));
%         evalStr = [];
%         for lTr = 1 : nTrs(mBase)
%             eval([cntrPeakNames{lTr} '=out{mBase, lTr}.cntrPeals;']);
%             %            eval([cntrPeakNames{lTr} '=out{mBase, lTr}.cntrPeaks;']);
%             evalStr = [evalStr, ',' cntrPeakNames{lTr}];
%         end
%         [tf{mBase}, locf{mBase}] = eval(['cellfun(DetectCmnPks' evalStr, ',', '''uniformoutput''' ', 0);']);
%         %        nStableCntrs{mBase} = 
%         STABLE_CNTRS{mBase}  = logical(cellfun(@sum, tf{mBase}));
%         for mTr = 1 : nTrs(mBase)
%             idx = tf{mBase}(STABLE_CNTRS{mBase});
%             cntrVertices = out{mBase, mTr}.cntrVertices(STABLE_CNTRS{mBase});
%             cntrPeaks = out{mBase, mTr}.cntrPeals(STABLE_CNTRS{mBase});
%             % commonCntrs.CntrPeaks{mBase, mTr} = out{mBase, mTr}.cntrPeaks(logical(STABLE_CNTRS{mBase}));
%             nStableCells = sum(STABLE_CNTRS{mBase});
%             nStableCntrs = cellfun(@sum, tf{mBase});
%             for kCell = 1 : nStableCells
%                 kIdx = idx{kCell};
%                 tempCntr = cntrVertices{kCell};
%                 cmnCntrs.cntrVertices{mBase, mTr} = tempCntr{kIdx};
%             end
%         end
%     end
    keyboard;
    f11 = offset;
    offset = offset(logical(sum(~cellfun(@isempty, offset), 2)), :);
    idx = nchoosek(1:size(offset, 2), 2);
    offset = reshape(offset(:, idx), [], 2);
    offset(logical(sum(cellfun(@isempty, offset), 2)), :) = [];
    pooledOffset = [];
    for ii = 1 : size(offset, 1)
        pooledOffset = [pooledOffset; offset{ii, 1}', offset{ii, 2}'];
    end
    xx = pooledOffset;
    xx(isnan(xx(:, 1)) | isnan(xx(:,2)), :) = [];
    fy = polyval([1, 0], xx(:, 1));
    resudials = vnorm([xx(:, 1), fy]' - xx');
    nPairs = size(pooledOffset, 1); 
    for kResample = 1 : nResample 
        % scramble cell ids
        pfr = [pooledOffset(randperm(nPairs), 1), pooledOffset(randperm(nPairs), 2)];
%         kRho = corrcoef(pfr, 'rows', 'pairwise');
%         rho(kResample) = kRho(1,2);
    end 
%     [cc, ii]=hist(rho, 1e3);    
%     bar(ii, cc, 'k');
%     hold on;
%     line([0.018, 0.018], ylim, 'color', 'r');
%     title(['p value = ' num2str(0.2442)], 'FontSize', 14);
    
    keyboard;

end

function funcHdl = GenDetectCmnPks(n, varargin)
% generates func hdl with variable no of input args for detecting common peaks across sessions
 
    [funcName, tolerence] = DefaultArgs(varargin, {'ismemberf', 2});
    str = [];
    for i = 1 : n    
        if i == 1, str = 'A1'; 
        else
            str = [str, ',A' num2str(i)];
        end
    end
    funcHdl = str2func(['@(' str ')' funcName '(' str ',' , '''row''', ',', '''tol''',  ',', num2str(tolerence), ')']);
end