function PoolOffset(varargin)

    [datasetType, roi, arena] = DefaultArgs(varargin, {'kenji', 'CA3', 'bigSquare'});
    switch datasetType
        case 'kenji'
          sK.roi = roi;
          sK.arena = arena;
          filebases = SearchKenji(sK);
          filebases = filebases(:, 1);
          
    end
    
    for kBase = 1 : length(filebases)
        trialNames = TrialNames(filebases{kBase}, datasetType, roi, arena);
        for mTr = 1 : length(trialNames)
            %            gt =  GenericTrial(filebases{kBase}, trialNames{mTr});
            try
                % load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(arena, roi), 'CCG.mat' ]);
                load(['~/data/analysis/kenji/', filebases{kBase}, '/', filebases{kBase}, '.', trialNames{mTr}, GenFiletag(arena, roi), 'CCG.mat' ]);
                offset{kBase, mTr} = fout.offset;
            catch err
            
            end     
kBase 
  end
    end
    f11 = offset;
    offset = offset(logical(sum(~cellfun(@isempty, offset), 2)), :);
    idx = nchoosek(1:size(offset, 2), 2);
    offset = reshape(offset(:, idx), [], 2);
    offset(logical(sum(cellfun(@isempty, offset), 2)), :) = [];
    pooledOffset = [];
    for ii = 1 : size(offset, 1)
        pooledOffset = [pooledOffset; offset{ii, 1}', offset{ii, 2}'];
    end

keyboard;
end
