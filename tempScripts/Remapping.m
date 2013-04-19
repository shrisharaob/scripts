function Remapping(filebase, varargin)
% Remapping(trial)
% 
    [arena, roi] = DefaultArgs(varargin, {[], {'bigSquare'}, {'CA3'}});
%     gt = GenericTrial();
    kenjiSearch.roi = roi;
    kenjiSearch.arena = arena;
    matches = SearchKenji(kenjiSearch);
    matches = matches(strcmp(matches(:, 1), filebase), :);
    nTrials = length(matches);
    refTr = 1;

    for kTr = 1 : nTrials
        gt = GenericTrial(filebase, matches{kTr, 2});
        if kTr == 1
            load([gt.paths.analysis, '.', gt.filebase, GenFiletag(roi, arena), '.commonClus.mat']);
        end
        fprintf(['\n trial :' gt.trialName ]);
        gt = gt.LoadPF;
        if ~isempty(gt.pfObject)
            %       roiClus = cell2mat(gt.GetRegionClu(roi));
            %roiPFUnits{kTr} = ismember(gt.pfObject.acceptedUnits, roiClus);
            %acceptedPFUnits{kTr} = gt.pfObject.acceptedUnits;
            roiPFPairs{kTr} = gt.pfObject.selectedPairs(ismember(gt.pfObject.selectedPairs, nchoosek(commonClus, 2), 'rows'), :);
            ratePk{kTr} = gt.pfObject.ratePk; % (ismember(gt.pfObject.acceptedUnits, trCluIdx{kTr}));
            pkDist{kTr} = gt.pfObject.pkDist;
        end
    end
   
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

    