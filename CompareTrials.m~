function CompareTrials(filebase, varargin)

    [arena, roi, cluIdx, IF_SmthRM, IF_REPORTFIG] = ...
        DefaultArgs(varargin, {{'bigSquare'}, {'CA1'}, [], 1, 0});
    kenjiSearch.roi = roi;
    kenjiSearch.arena = arena;
    matches = SearchKenji(kenjiSearch);
    matches = matches(strcmp(matches(:, 1), filebase), :);
    nTrials = size(matches, 1);
    refTr = 1;
    commonClus = cluIdx;
    if nTrials == 1, return; end;
    for kTr = 1 : nTrials
        gt = GenericTrial(filebase, matches{kTr, 2});
        if kTr == 1 & isempty(cluIdx)
            load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
        end 
        if length(commonClus) == 1, return; end; % if only one common unit across the trials
        fprintf(['\n trial :' gt.trialName ]);
        gta{kTr} = gt.LoadPF;
        sparsity(:, kTr) = gta{kTr}.pfObject.sparsity(ismember(gta{kTr}.pfObject.acceptedUnits, commonClus));
        if kTr == refTr
            [~, sRank] = sort(sparsity(:, refTr), 'ascend');
            commonClus = commonClus(sRank);
        end
    end
    sparsity = sparsity(sRank, :); 
    hFig = figure;
    for lClu = 1 : length(commonClus)
        if IF_SmthRM
            for mTr = 1 : nTrials
                subplot(1, nTrials, mTr);
                curgt = gta{mTr};
                curgt.pfObject.PlotRateMaps(1, 0, 1, [],[],[], commonClus(lClu));
                axis square;
                srm = gta{mTr}.pfObject.smoothRateMap(:,:,ismember(gta{mTr}.pfObject.acceptedUnits, commonClus(lClu)));
                mapEntropy(lClu, mTr) = Entropy(srm);
                
            end
        end
        
        if IF_REPORTFIG
            reportfig(hFig, [mfilename, filebase, GenFiletag(roi,arena)], 0, [filebase, '    sparsity : ' num2str(sparsity(lClu, :)), 'entr  ' num2str(mapEntropy(lClu, :))]);
        end
        clf;
    end
end

