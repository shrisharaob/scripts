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
    refTr = 1;
    for kTr = 1 : nTrials
        if isempty(filebase)
            gt = GenericTrial(matches{kTr, 1}, matches{kTr, 2});
        else 
            gt = GenericTrial(filebase, matches{kTr, 2});
        end
        fprintf(['\n trial :' gt.trialName ]);
        gt = gt.LoadPF;
       
        if ~isempty(gt.pfObject)
            roiClus = cell2mat(gt.GetRegionClu(roi));
            roiPFUnits{kTr} = ismember(gt.pfObject.acceptedUnits, roiClus);
            acceptedPFUnits{kTr} = gt.pfObject.acceptedUnits;
            roiPFPairs{kTr} = gt.pfObject.selectedPairs(ismember(gt.pfObject.selectedPairs, nchoosek(roiClus, 2), 'rows'), :);
            ratePk{kTr} = gt.pfObject.ratePk; % (ismember(gt.pfObject.acceptedUnits, trCluIdx{kTr}));
            pkDist{kTr} = gt.pfObject.pkDist;
        end
    end
%     trCluIdx{kTr}= clus(ismember(clus, linCluIdx(gt.pfObject.idealPFUnits)));

    for kTr = setdiff(1 : nTrials, refTr)
        if ~isempty(roiPFUnits{kTr})
            refPk = ratePk{refTr};
            % is a cell which is in roi & is a pyr & meets pf criterion
            refPFUnitIdx = roiPFUnits{refTr} & ismember(acceptedPFUnits{refTr}, acceptedPFUnits{kTr});
            refPk = refPk(refPFUnitIdx);
            kPk = ratePk{kTr};
            kPFUnitIdx = roiPFUnits{kTr} & ismember(acceptedPFUnits{kTr}, acceptedPFUnits{refTr});
            kPk = kPk(kPFUnitIdx);
            refPkDist = pkDist{refTr};
            refPkDist = refPkDist(ismember(roiPFPairs{refTr}, roiPFPairs{kTr}, 'rows'));
            kPkDist = pkDist{kTr};
            kPkDist = kPkDist(ismember(roiPFPairs{kTr}, roiPFPairs{refTr}, 'rows'));
            figure(kTr)
%             clf;
            subplot(1, 2, 1);
            if length(kPk) >= 2
                hist2([MakeUniformDistr(refPk), MakeUniformDistr(kPk)]);
            end
            subplot(1, 2, 2)
            if length(kPkDist) >= 2
                hist2([MakeUniformDistr(refPkDist), MakeUniformDistr(kPkDist)]);
            end
%             xlabel(['trial :', ]);
%             ylabel(['trial# 
%             waitforbuttonpress;
%             clf;
        end
    end
end

    