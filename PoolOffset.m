function PoolOffset(varargin)

    [datasetType, roi, arena, nResample] = DefaultArgs(varargin, {'kenji', 'CA3', 'bigSquare', 1e2});
    switch datasetType
        case 'kenji'
          sK.roi = roi;
          sK.arena = arena;
          filebases = SearchKenji(sK);
          filebases = filebases(:, 1);
          
    end
    
    for kBase = 1 : length(filebases)
        trialNames = TrialNames(filebases{kBase}, datasetType, roi, arena);
        fprintf([filebases{kBase} '\n']);
        for mTr = 1 : length(trialNames)
            fprintf(['\t' trialNames{mTr} '\n']);
            gt =  GenericTrial(filebases{kBase}, trialNames{mTr});
            try
                % load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(arena, roi), 'CCG.mat' ]);
                load(['~/data/analysis/kenji/', filebases{kBase}, '/', filebases{kBase}, '.', trialNames{mTr}, GenFiletag(arena, roi), 'CCG.mat' ]);
                offset{kBase, mTr} = abs(fout.offset);
                commonClus = gt.LoadCommonClus;

                out{kBase, mTr} = gt.MultiPeakPFDistance(nchoosek(commonClus, 2));
            catch err
                
            end     
            %kBase 
        end
    end
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
