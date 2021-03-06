function catPeriods = CatTrials(filebase, roi, arena, state, varargin )
% out = CatTrials(filebase, roi, arena)
    
    [outFs] = DefaultArgs(varargin, {0});
    out  = [];
    [matches, fbRoi] =  SearchKenji(filebase);
    if ~any(cellfun(@strcmp, fbRoi, repmat({roi}, size(fbRoi)))); return; end
    trialNames = matches(cellfun(@strcmp, matches(:, 3), repmat({arena}, size(matches(:, 3)))), 2);
    trialPeriods = [];
    for mTr = 1 : length(trialNames)
        gt = GenericTrial(filebase, trialNames{mTr});
        trialPeriods = [trialPeriods; gt.LoadStatePeriods(state, 0, 0)];
    end
    statePeriods = load([gt.paths.data, gt.filebase '.sts.', state]);
    catPeriods = IntersectRanges(trialPeriods, statePeriods);
    if outFs ~= 0
        catPeriods = ConvertFs(catPeriods, gt.lfpSampleRate, outFs);
    end
end