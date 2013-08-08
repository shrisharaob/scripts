function out = Remapping(filebase, datasetType, varargin)
% Remapping(filebase, arena, roi, cluIdx, IF_PLOT)

    [roi, arena, cluIdx, IF_PLOT, IF_REPORTFIG, sparsityThresh] = ...
        DefaultArgs(varargin, {{'CA3'}, {'bigSquare'}, [], 0, 0, 0.1});
    out = {};
%   kenjiSearch.roi = roi;
%   kenjiSearch.arena = arena;
%   matches = SearchKenji(kenjiSearch);
%   matches = matches(strcmp(matches(:, 1), filebase), :);
    trialNames = TrialNames(filebase, datasetType, roi, arena);
    nTrials = length(trialNames);
    if nTrials < 2, return; end
    nTrPairs = nchoosek(nTrials, 2);
    refTr = 1;
    for kTr = 1 : nTrials
        gt = GenericTrial(filebase, trialNames{kTr});
        gt.LoadPF;
        cluIdx = gt.pfObject.acceptedUnits(gt.pfObject.sparsity < sparsityThresh);
        if kTr == 1 
            %load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
            commonClus = FindCommonClu(gt.filebase, [], roi, arena, sparsityThresh);
            commonClus = intersect(commonClus, cluIdx);
        end
        if length(commonClus) < 2, return; end
        fprintf(['\n trial :' gt.trialName ]);
        if ~isempty(gt.pfObject)
            roiClus = cell2mat(gt.GetRegionClu(roi));
            roiPFUnits{kTr} = ismember(gt.pfObject.acceptedUnits, roiClus);
            acceptedPFUnits{kTr} = gt.pfObject.acceptedUnits;
            if length(acceptedPFUnits{kTr}) > 1, tempPairs = nchoosek(acceptedPFUnits{kTr}, 2); else, return; end
            %roiPFPairs{kTr} = gt.pfObject.selectedPairs(ismember(gt.pfObject.selectedPairs, nchoosek(commonClus, 2), 'rows'), :);
            roiPFPairs{kTr} = tempPairs(ismember(tempPairs, nchoosek(commonClus, 2), 'rows'), :);
            ratePk{kTr} = gt.pfObject.ratePk; % (ismember(gt.pfObject.acceptedUnits, trCluIdx{kTr}));
            pkDist{kTr} = gt.pfObject.pkDist;
        end
    end
    out.pkRateMat = [];
    out.pkDistMat = [];
    for kTr = setdiff(1 : nTrials, refTr)
        if ~isempty(roiPFUnits{kTr}) & ~isempty(roiPFPairs{kTr})
            refPk = ratePk{refTr};
            % isempty a cell which is in roi & is a pyr & meets pf criterion
            refPFUnitIdx = roiPFUnits{refTr} & ismember(acceptedPFUnits{refTr}, acceptedPFUnits{kTr});
            refPk = refPk(refPFUnitIdx);
            kPk = ratePk{kTr};
            kPFUnitIdx = roiPFUnits{kTr} & ismember(acceptedPFUnits{kTr}, acceptedPFUnits{refTr});
            kPk = kPk(kPFUnitIdx);
            refPkDist = pkDist{refTr};
            refPkDist = refPkDist(ismember(roiPFPairs{refTr}, roiPFPairs{kTr}, 'rows'));
            kPkDist = pkDist{kTr};
            kPkDist = kPkDist(ismember(roiPFPairs{kTr}, roiPFPairs{refTr}, 'rows'));
            %            out.pkDist{kTr} = [refPkDist, kPkDist];
            out.pkRateMat = [out.pkRateMat; refPk', kPk'];
            out.pkDistMat = [out.pkDistMat; refPkDist', kPkDist'];

            if IF_PLOT
                figure(kTr)
                subplot(nTrPairs, 2, 1);
                if length(kPk) >= 2
                    %                    hist2([MakeUniformDistr(refPk), MakeUniformDistr(kPk)]);
                    plot(refPk, kPk, '+');
                    title('rate peaks');
                    xlabel(trialNames{refTr});
                    ylabel(trialNames{kTr});
                    axis square;
                end
                subplot(nTrPairs, 2, 2)
                if length(kPkDist) >= 2 & ~isempty(kPk)
                    %hist2([MakeUniformDistr(refPkDist), MakeUniformDistr(kPkDist)]);
                    plot(refPk, kPkDist, '+');
                    title('peak distances');d
                    xlabel(trialNames{refTr});
                    ylabel(trialNames{kTr});
                    axis square;
                end
                if IF_REPORTFIG
                    reportfig(kTr, [mfilename, GenFiletag(roi, arena), datasetType], 0, []);
                end
            end
        end
    end
end

    
% plot(out.poolArray(:,1), out.poolArray(:,2),'+k')
%  hold  on
%  xlabel('Rate Peak Trial A')
%  xlim([0, 40])
%  ylabel('Rate Peak Trial B ')
%  ylim([0, 40])
%  line(xlim, ylim, 'color', 'k')
% [r,p] = corrcoef(out.poolArray)
%  text(6, 33, ['r = ', num2str(r(1, 2), '%.2f')], 'FontSize', 12, 'FontWeight', 'b')
%  ProcessFigure(1, '~/thesis/figures/report/RateRemapping_CA3_bigSquare')