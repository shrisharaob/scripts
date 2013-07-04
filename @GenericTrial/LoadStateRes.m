function [res, clu, varargout] = LoadStateRes(gt, varargin)
% [res, clu, pos, posStatePeriods] = LoadStateRes(gt, varargin)
% [state, IF_INGOODPOS, fs, posInPeriods] = {0, 0}
% IF_INGOODPOS : logical, res only in good pos periods loaded
% fs : convert res to fs
% loads clu res for the specified state

    [state, IF_INGOODPOS, fs, posInPeriods, IF_SQUASH] = DefaultArgs(varargin, {'RUN', 1, 0, [], 0});

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
        if ~strcmp(state, 'trajEvnts'),
            statePeriods = load([gt.paths.data, gt.filebase '.sts.', state]); % @lfp fs
        end
        switch state
          case 'SWS'
             [res, resIdx] = SelectPeriods(gt.res, ConvertFs(statePeriods, gt.lfpSampleRate, gt.sampleRate), 'd', 1, IF_SQUASH);
             clu = gt.clu(resIdx);
             varargout = {[], []};
          case 'trajEvnts'
            sts = gt.TrajectoryEvents(1);
            [res, resIdx] = SelectPeriods(gt.res, sts, 'd', 1, 1);
            clu = gt.clu(resIdx);
            varargout = {[], []};
          case 'RUN'
            markerNo = 1;
            statePeriods = IntersectRanges(statePeriods, gt.trialPeriods);
            trialStartTime_pos =  ConvertFs(gt.trialPeriods(1, 1), gt.lfpSampleRate, gt.trackingSampleRate);
            if IF_INGOODPOS
                posStatePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.trackingSampleRate); %round(statePeriods .* gt.trackingSampleRate  ./ gt.lfpSampleRate) + 1;
                                                                                                    % statePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.sampleRate);
                posStatePeriods = IntersectRanges(posStatePeriods, gt.goodPosPeriods + trialStartTime_pos); 
                statePeriods = ConvertFs(posStatePeriods, gt.trackingSampleRate, gt.sampleRate);
            else
                statePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.sampleRate);
                posStatePeriods = ConvertFs(statePeriods, gt.sampleRate, gt.trackingSampleRate);
            end


            if isempty(posInPeriods) 
                [res, resInd] = SelectPeriods(gt.res, statePeriods, 'd', 1, IF_SQUASH);
                clu = gt.clu(resInd);
                if fs > 0
                    res = ConvertFs(res, gt.sampleRate, fs);
                end
                varargout{1} = SelectPeriods(sq(gt.position(:,markerNo,:)), posStatePeriods - trialStartTime_pos, 'c');
            else
                for kInterval = 1 : size(posInPeriods, 1)
                    posIntervals = RecenterPeriods(posInPeriods);
                    pos{kInterval} = SelectPeriods(sq(gt.position(:, markerNo, :)), IntersectRanges(posStatePeriods - trialStartTime_pos, posIntervals(kInterval, :)), 'c');
                    kRes = ConvertFs(gt.res, gt.sampleRate, gt.trackingSampleRate);
                    [res{kInterval}, kResIndx] = SelectPeriods(kRes, IntersectRanges(posStatePeriods, posInPeriods(kInterval, :)), 'd', 1, 1);
                    clu{kInterval} = gt.clu(kResIndx);
                    if fs > 0 && fs ~= gt.trackingSampleRate
                        res{kInterval} = ConvertFs(res{kInterval}, gt.trackingSampleRate, fs);
                    end
                end
                varargout{1} = pos;
            end
            varargout{2} = posStatePeriods;
        end
    end
end