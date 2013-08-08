function out = CommonCluStats(gt, roi, arena, varargin)

    sparsityThresh = DefaultArgs(varargin, {.1});
    gt.LoadPF;
    out = [];
    matches = SearchKenji(gt.filebase);
    arenaNames = matches(:, 3);
    if ~iscell(arena), arena = cellstr(arena); end
    if length(arena) == 1
        normFactor = sum(cell2mat(cellfun(@strcmp, arenaNames, repmat({char(arena)}, length(arenaNames), 1), 'UniformOutput', 0)));
    else
        normFactor = sum(cell2mat(cellfun(@strcmp, arenaNames, repmat({char(arena{1})}, length(arenaNames), 1), 'UniformOutput', 0)) | cell2mat(cellfun(@strcmp, arenaNames, repmat({char(arena{2})}, length(arenaNames), 1), 'UniformOutput', 0)));
    end
    try 
        if length(arena) == 2
            curArena = SearchKenji(gt.trialName);
            if strcmp(curArena{2}, 'linear'), sparsityThresh = .35;  end
            nCommonClu = length(intersect(gt.LoadCommonClus(roi, arena), gt.pfObject.acceptedUnits(gt.pfObject.sparsity <= sparsityThresh)));
        else 
            switch char(arena)
              case 'bigSquare'
                mpd = gt.MultiPeakPFDistance(roi, arena{1});
                outTemp = histc(cellfun(@length, mpd.cntrVertices), 1:10);
                nCommonClu = outTemp(1);
              case 'linear'
                nCommonClu = length(intersect(gt.LoadCommonClus(roi, arena), gt.pfObject.acceptedUnits(gt.pfObject.sparsity <= 0.3)));
            end
        end
        %out.commonclu = gt.LoadCommonClus(roi, arena);
        out.nCmnClus = length(nCommonClu) ./ normFactor;
    catch err
        fprintf([err.message '\n']);
    end
end

