function PoolOffset(varargin)

    [type, datasetType, roi, arena, nResample] = ...
        DefaultArgs(varargin, {'load', 'kenji', 'CA3', {'bigSquare'}, 1e2});
    
    switch  type
      case 'compute'
        switch datasetType
          case 'kenji'
            sK.roi = roi;
            sK.arena = arena;
            filebases = SearchKenji(sK);
            filebases = unique(filebases(:, 1));
        end
        trNo = 0;
        for lArena = 1 : length(arena)
            for kBase = 1 : length(filebases)
                trialNames = TrialNames(filebases{kBase}, datasetType, roi, arena{lArena});
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
                            tcntrs{kBase, mTr} = []; %struct('cntrVertices', [], 'cntrPeaks', [], 'cluId', []);
                        elseif length(commonClus) < 2
                            tcntrs{kBase, mTr} = []; %struct('cntrVertices', [], 'cntrPeaks', [], 'cluId', []);
                        else
                            tcntrs{kBase, mTr} = gt.MultiPeakPFDistance( roi, arena, nchoosek(commonClus, 2));
                            if isempty({tcntrs{kBase, mTr}.cluId}), tcntrs{kBase, mTr} = []; end
                        end
                        trNo = trNo + 1;
                        trNames{trNo} = trialNames{mTr}; 
                    catch err
                        fprintf(['\n error ' trialNames{mTr}, '\n']);
                    end     
                end
            end
            cntrs{lArena} = tcntrs;
            clear tcntrs;
        end
        save(['~/data/analysis/', datasetType, '/Contours', GenFiletag(roi, arena), 'mat'], 'cntrs', 'filebases');

        %%
      case 'load'
        load(['~/data/analysis/', datasetType, '/Contours', GenFiletag(roi, arena), 'mat'], 'cntrs', 'filebases');
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
        end 
        keyboard;
    end
end
