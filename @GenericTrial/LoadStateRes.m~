function [res, clu] = LoadStateRes(gt, varargin)
%[clu, res] = LoadStateRes(gt, varargin)
% [IF_INGOODPOS, fs] = {0, 0}
% IF_INGOODPOS : logical, res only in good pos periods loaded
% fs : convert res to fs
% loads clu res for the specified state

    [state, IF_INGOODPOS, fs] = DefaultArgs(varargin, {'RUN', 0, 0});

    if isempty(gt.res), gt = gt.LoadCR; end

    switch gt.datasetType
      case  'MTA'
        mtaTrial = MTATrial(gt.filebase, [], gt.trialName);
        if strcmp(state, 'RUN')
            state = 'walk'; 
        else
            error('\n no %s state in filebase %s', state, gt.filebase); 
        end
        statePeriods = mtaTrial.statePeriods(state); %% state periods starts with 0
        markerNo = 7;
        %posStatePeriods = round(statePeriods .* gt.trackingSampleRate ./ gt.lfpSampleRate) + 1;
        statePeriods = round(statePeriods .* gt.sampleRate ./ gt.lfpSampleRate) + 1;
        trialStartTime_pos = 0;
      case 'default'
        statePeriods = load([gt.paths.data, gt.filebase '.sts.', state]);
        % posStatePeriods = round(statePeriods .* gt.trackingSampleRate  ./ gt.lfpSampleRate) + 1;
        statePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.sampleRate);
        trialStartTime_pos = 0;
        markerNo = 1;
      case  'kenji'
        statePeriods = load([gt.paths.data, gt.filebase '.sts.', state]); % @lfp fs
        statePeriods = IntersectRanges(statePeriods, gt.trialPeriods);
       if IF_INGOODPOS
           posStatePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.trackingSampleRate); %round(statePeriods .* gt.trackingSampleRate  ./ gt.lfpSampleRate) + 1;
           % statePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.sampleRate);
           trialStartTime_pos =  ConvertFs(gt.trialPeriods(1, 1), gt.lfpSampleRate, gt.trackingSampleRate);
           posStatePeriods = IntersectRanges(posStatePeriods, gt.goodPosPeriods + trialStartTime_pos); 
           statePeriods = ConvertFs(posStatePeriods, gt.trackingSampleRate, gt.sampleRate);
       else
           statePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.sampleRate);
       end
    end

    [res, resInd] = SelectPeriods(gt.res, statePeriods, 'd', 1, 1);
    clu = gt.clu(resInd);
    if fs > 0
        res = ConvertFs(res, gt.sampleRate, fs);
    end
end
