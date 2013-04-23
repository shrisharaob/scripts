function out = PairRemapping(filebase, varargin)
% PairRemapping(filebase, cluIdx, arena, roi, IF_PLOT, IF_REPORTFIG)
%
    [cluIdx, arena, roi, IF_PLOT, IF_REPORTFIG] = ...
        DefaultArgs(varargin, {[], {'bigSquare'}, {'CA3'}, 0, 0});
    out = {};
    kenjiSearch.roi = roi;
    kenjiSearch.arena = arena;
    matches = SearchKenji(kenjiSearch);
    matches = matches(strcmp(matches(:, 1), filebase), :);
    nTrials = size(matches, 1);
    refTr = 1;
    commonClus = cluIdx;
    %    if nTrials == 1, return; end;
    for kTr = 1 : nTrials
        gt = GenericTrial(filebase, matches{kTr, 2});
        if kTr == 1  & isempty(cluIdx)
            load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
        end
        if length(commonClus) == 1, return; end; % if only one common unit across the trials
        fprintf(['\n trial :' gt.trialName ]);
        gt = gt.LoadPF;
        if ~isempty(gt.pfObject)
            %            roiPFPairs{kTr} = gt.pfObject.accepted(ismember(gt.pfObject.selectedPairs, nchoosek(commonClus, 2), 'rows'), :);
            idealCmnPFUnits = gt.pfObject.acceptedUnits(ismember(gt.pfObject.acceptedUnits, gt.pyrCluIdx(gt.pfObject.idealPFUnits)));
            ratePk{kTr} =  gt.pfObject.ratePk(ismember(gt.pfObject.acceptedUnits, commonClus)); %idealCmnPFUnits(ismember(idealCmnPFUnits, commonClus))));
            pkDist{kTr} = gt.pfObject.pkDist(ismember(gt.pfObject.acceptedUnits, commonClus));  %idealCmnPFUnits(ismember(idealCmnPFUnits, commonClus))));
        end
    end
    out.pkDist = [];
    if nTrials == 1
        out.pkDist = pkDist{1};
    end
    for kTr = setdiff(1 : nTrials, refTr)
        if ~isempty(commonClus)
            refPk = ratePk{refTr};
            kPk = ratePk{kTr};
            refPkDist = pkDist{refTr};
            kPkDist = pkDist{kTr};
            out.pkDist = [out.pkDist; refPkDist', kPkDist'];
            if IF_PLOT
                hFig = figure(kTr);
                subplot(1, 2, 1);
                plot(refPk, kPk, '*k');
                grid on;
                hold on;
                axis square;
                line([0, max(max(xlim), max(ylim))], [0, max(max(xlim), max(ylim))], 'color', 'k');
                title('peak rate');
                xlabel(['reference trial :' matches{refTr, 2}]);
                ylabel(['trial id:' matches{kTr, 2}]);
                subplot(1, 2, 2)
                plot(refPkDist, kPkDist, '*k');
                hold on;
                grid on;
                line([0, max(max(xlim), max(ylim))], [0, max(max(xlim), max(ylim))], 'color', 'k');
                axis square;
                title('euclidean dist between peaks');
                xlabel(['reference trial :' matches{refTr, 2}]);
                ylabel(['trial id:' matches{kTr, 2}]);
                drawnow;
                if IF_REPORTFIG
                    reportfig(hFig, [mfilename, GenFiletag(roi, arena) '.rates_&_pkDist'], 0, filebase);
                end
            end
        end
    end
end


