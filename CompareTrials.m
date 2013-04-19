function CompareTrials(filebase, varargin)

    [arena, roi, IF_SmthRM, IF_REPORTFIG] = DefaultArgs(varargin, {{'bigSquare'}, {'CA3'}, 1, 1});
    kenjiSearch.roi = roi;
    kenjiSearch.arena = arena;
    matches = SearchKenji(kenjiSearch);
    matches = matches(strcmp(matches(:, 1), filebase), :);
    nTrials = size(matches, 1);
    refTr = 1;
    if nTrials == 1, return; end;
    for kTr = 1 : nTrials
        gt = GenericTrial(filebase, matches{kTr, 2});
        if kTr == 1
            load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
        end 
        if length(commonClus) == 1, return; end; % if only one common unit across the trials
        fprintf(['\n trial :' gt.trialName ]);
        gta{kTr} = gt.LoadPF;
    end
    hFig = figure;
    for lClu = 1 : length(commonClus)
        if IF_SmthRM
            for mTr = 1 : nTrials
                subplot(1, nTrials, mTr);
                curgt = gta{mTr};
                curgt.pfObject.PlotRateMaps(0, 0, 1, [],[],[], commonClus(lClu));
                axis square;
            end
        end
%         waitforbuttonpress;
    if IF_REPORTFIG
        reportfig(hFig, [mfilename, filebase, GenFiletag(roi, arena)], 0, filebase);
    end
    end
end