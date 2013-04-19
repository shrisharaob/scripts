function PairRemapping(filebase, varargin)
    % Remapping(trial)
    %
    [arena, roi, IF_REPORTFIG] = DefaultArgs(varargin, {{'bigSquare', 'linear'}, {'CA3'}, 1});
    %     gt = GenericTrial();
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
        gt = gt.LoadPF;
        if ~isempty(gt.pfObject)
            %       roiClus = cell2mat(gt.GetRegionClu(roi));
            %roiPFUnits{kTr} = ismember(gt.pfObject.acceptedUnits, roiClus);
            %acceptedPFUnits{kTr} = gt.pfObject.acceptedUnits;
            roiPFPairs{kTr} = gt.pfObject.selectedPairs(ismember(gt.pfObject.selectedPairs, nchoosek(commonClus, 2), 'rows'), :);
            ratePk{kTr} =  gt.pfObject.ratePk(ismember(gt.pfObject.acceptedUnits, commonClus));
            pkDist{kTr} = gt.pfObject.pkDist(ismember(gt.pfObject.acceptedUnits, commonClus));
        end
    end
    for kTr = setdiff(1 : nTrials, refTr)
        if ~isempty(commonClus)
            refPk = ratePk{refTr};
            % is a cell which is in roi & is a pyr & meets pf criterion
            %            refPFUnitIdx = roiPFUnits{refTr} & ismember(acceptedPFUnits{refTr}, acceptedPFUnits{kTr});
            kPk = ratePk{kTr};
            refPkDist = pkDist{refTr};
            kPkDist = pkDist{kTr};
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
            %             if length(kPk) >= 2
            %                 hist2([MakeUniformDistr(refPk), MakeUniformDistr(kPk)]);
            %             end
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
        
        %             if length(kPkDist) >= 2
        %                 hist2([MakeUniformDistr(refPkDist), MakeUniformDistr(kPkDist)]);
        %             end
        %             xlabel(['trial :', ]);
        %             ylabel(['trial#
        %             waitforbuttonpress;
        %             clf;
        
    end
end


