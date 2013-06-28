function statePeriods = LoadStatePeriods(gt, varargin)
% statePeriods = LoadStatePeriods(gt, varargin)
% [state, fs, IF_INGOODPOS]
% loads the specified state periods @lfp fs

    [state, fs, IF_INGOODPOS] = DefaultArgs(varargin, {'RUN', 0, 1});

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
        posStatePeriods = round(statePeriods .* gt.trackingSampleRate ./ gt.lfpSampleRate) + 1;
        statePeriods = round(statePeriods .* gt.sampleRate ./ gt.lfpSampleRate) + 1;
        trialStartTime_pos = 0;
      case 'default'
        statePeriods = load([gt.paths.data, gt.filebase '.sts.', state]);
        posStatePeriods = round(statePeriods .* gt.trackingSampleRate  ./ gt.lfpSampleRate) + 1;
        statePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.sampleRate);
        trialStartTime_pos = 0;
        markerNo = 1;
      case  'kenji'
        statePeriods = load([gt.paths.data, gt.filebase '.sts.', state]); % @lfp fs
        switch state
          case 'RUN'
            markerNo = 1;
            statePeriods = IntersectRanges(statePeriods, gt.trialPeriods);
          case 'traj'
            statePeriods = gt.TrajectoryEvents(state);
            return;
          case 'SWS'
            if fs ~= 0
                statePeriods = ConvertFs(statePeriods, gt.trackingSampleRate, fs);
            end
            return;
        end
    end
    
    if IF_INGOODPOS
        posStatePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.trackingSampleRate); %round(statePeriods .* gt.trackingSampleRate  ./ gt.lfpSampleRate) + 1;
        trialStartTime_pos =  ConvertFs(gt.trialPeriods(1, 1), gt.lfpSampleRate, gt.trackingSampleRate);
        posStatePeriods = IntersectRanges(posStatePeriods, gt.goodPosPeriods + trialStartTime_pos); 
    end
    
    if fs ~= 0
        statePeriods = ConvertFs(posStatePeriods, gt.trackingSampleRate, fs);
    elseif fs == 0 & IF_INGOODPOS
        statePeriods = ConvertFs(posStatePeriods, gt.trackingSampleRate, gt.lfpSampleRate);
    end
end




