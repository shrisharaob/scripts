function SequenceMem(gt, varargin)
% function SequenceMem(gt, varargin)
% [arena, roi, IF_REPORTFIG, type] 
% {'bigSquare'}, {'CA3'}, 1, 'display'
% script for decoding position of the rat using the rate maps from
% other trials 
% the decoder stisfies strict Cramer-Rao bound under conditions of
% poission spikes stats and independence of cells
% Zhang 1998, J. Neuro phys

    [ arena, roi, IF_REPORTFIG, type] = ...
        DefaultArgs(varargin, { {'bigSquare'}, {'CA3'}, 1, 'display'});
  
    filebase = gt.filebase;
    trialNames = TrialNames(gt.filebase, gt.datasetType, roi, arena);
    nTrials = length(trialNames);
    refTr = nTrials;
    switch gt.datasetType
        case 'MTA'
          roi = 'CA1';
          arena = 'cof';
    end
    switch type
      case 'compute'
        if nTrials == 1, return; end;
        for kTr = 1 : nTrials
            gt = GenericTrial(filebase, trialNames{kTr});
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
                save([gt.paths.analysis, filebase,  GenFiletag(roi, arena), mfilename, trialNames{mTr}], 'pos', 'predErr');
                drawnow;
                if IF_REPORTFIG
                    reportfig(gcf, [mfilename, GenFiletag(roi, arena) 'rates_&_pkDist'], 0, [filebase ':: refTrial - ' char(SearchKenji(trialNames{refTr})), '(' trialNames{refTr} ')  trial - ' ...
                                        char(SearchKenji(trialNames{refTr})) '(' trialNames{mTr} ')'] );
                    close gcf;
                end
            end
        end
      case 'display'
        gt = GenericTrial(filebase, trialNames{refTr});
        for mTr = 1 : nTrials
            %            load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
            load([gt.paths.analysis, filebase,  GenFiletag(roi, arena), mfilename, trialNames{mTr}], 'pos', 'predErr');
            h = figure;
            plot(predErr);
            title([num2str(binSize * 1e3), 'ms bins']);
            if IF_REPORTFIG
                reportfig(gcf, [mfilename, GenFiletag(roi, arena) 'rates_&_pkDist'], 0, [filebase ':: refTrial - ' char(SearchKenji(trialNames{refTr})), '(' trialNames{refTr} ')  trial - ' ...
                                    char(SearchKenji(trialNames{refTr})) '(' trialNames{mTr} ')'] );
                close(h);
            end
        end
    end
    close all;
end
