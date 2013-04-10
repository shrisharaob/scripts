function Remapping(varargin)
% Remapping(trial)
% 
    [filebase, arena, roi] = DefaultArgs(varargin, {[], {'bigSquare'}, {'CA3'}});
%     gt = GenericTrial();
    kenjiSearch.roi = roi;
    kenjiSearch.arena = arena;
    matches = SearchKenji(kenjiSearch);
    if ~isempty(filebase)
%         gt = genericTfilebase;
        matches = matches(strcmp(matches(:, 1), filebase), :);
    end
    nTrials = length(matches);
    for kTr = 1 : nTrials
        if isempty(filebase)
            gt = GenericTrial(matches{kTr, 1}, matches{kTr, 2});
        else 
            gt = GenericTrial(filebase, matches{kTr, 2});
        end
        gt = gt.LoadPF;
        if ~isempty(gt.pfObject)
            ratePk{kTr} = gt.pfObject.ratePk;
            clus = cell2mat(gt.GetRegionClu(roi));
            trCluIdx{kTr}= clus(ismember(clus, gt.pfObject.idealPFUnits));
        end
    end
    refTr = 1;
    for kTr = setdiff(1 : nTrials, refTr)
        if ~isempty(trCluIdx{kTr})
            hist2(MakeUniformDistr(ratePk{refTr}), MakeUniformDistr(ratePk{kT}));
            waitforbuttonpress;
            clf;
        end
    end
end

    