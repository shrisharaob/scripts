function SequenceMem(filebase, varargin)
% script for decoding position of the rat using the rate maps from
% other trials 
% the decoder stisfies strict Cramer-Rao bound under conditions of
% poission spikes stats and independence of cells

[arena, roi, IF_REPORTFIG, type] = DefaultArgs(varargin, {{'bigSquare'}, {'CA3'}, 1, 'display'});

kenjiSearch.roi = roi;
kenjiSearch.arena = arena;
matches = SearchKenji(kenjiSearch);
matches = matches(strcmp(matches(:, 1), filebase), :);
nTrials = size(matches, 1);
refTr = nTrials;
switch type
  case 'compute'
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
            smRatemaps{kTr} = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus));
        end
        if kTr == refTr
            refGt = gt;
            if isempty(refGt.res)
                refGt = refGt.LoadCR;
            end
        end
    end
    for mTr = 1 : nTrials
        if ~isempty(commonClus)
            [pos, predErr] = DecodePos(refGt, smRatemaps{mTr}, commonClus, 200e-3);
            save([gt.paths.analysis, filebase,  GenFiletag(roi, arena), mfilename, matches{mTr, 2}], 'pos', 'predErr');
            drawnow;
            if IF_REPORTFIG
                reportfig(gcf, [mfilename, GenFiletag(roi, arena) 'rates_&_pkDist'], 0, [filebase ':: refTrial - ' char(SearchKenji(matches{refTr, 2})), '(' matches{refTr, 2} ')  trial - ' ...
                                    char(SearchKenji(matches{refTr, 2})) '(' matches{mTr, 2} ')'] );
                close gcf;
            end
        end
    end
  case 'display'
   gt = GenericTrial(filebase, matches{refTr, 2});
  for mTr = 1 : nTrials
      %            load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
            load([gt.paths.analysis, filebase,  GenFiletag(roi, arena), mfilename, matches{mTr, 2}], 'pos', 'predErr');
            h = figure;
            plot(predErr);
            title([num2str(binSize * 1e3), 'ms bins']);
            if IF_REPORTFIG
                reportfig(gcf, [mfilename, GenFiletag(roi, arena) 'rates_&_pkDist'], 0, [filebase ':: refTrial - ' char(SearchKenji(matches{refTr, 2})), '(' matches{refTr, 2} ')  trial - ' ...
                                    char(SearchKenji(matches{refTr, 2})) '(' matches{mTr, 2} ')'] );
                close(h);
            end
  end
end
close all;
end
