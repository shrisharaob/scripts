function out = ComparePVs(filebase, varargin)

    [roi, arena] = DefaultArgs(varargin, {{'CA3'}, {'bigSquare'}});
    filetag = [];
    for mRoi = 1 : length(roi)
        if mRoi == 1
            filetag = ['.', char(roi{1})];
        end
        if mRoi > 1
            filetag = [filetag, '.', char(roi{mRoi})];
        end
    end
    for lArena = 1 : length(arena)
        filetag = [filetag, '.', char(arena{lArena})];
    end
    sK.roi = roi;
    sK.arena = arena;
    matches = SearchKenji(sK);
    matches = matches(strcmp(matches(:, 1), filebase), :);
    trialNames = matches(:, 2);
    pvNames = genvarname(cellstr(repmat('pv', length(trialNames), 1)));
    rvNames = genvarname(cellstr(repmat('rv', length(trialNames), 1)));
    %    distNames = genVarNames(cellstr(repmat('rv', length(trialNames), 1)));
    for mTr = 1 : length(trialNames)
        load(['~/data/analysis/kenji/' filebase, '/', filebase, '.', trialNames{mTr}, filetag, 'PopVecTimeCourse.mat']);
        eval([pvNames{mTr} '= popVec;']);
        eval([rvNames{mTr} '= refVector;']);
    end
keyboard;
    %    distance = norm(rv - rv2);    
end
