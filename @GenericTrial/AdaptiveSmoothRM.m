function smoothRateMap = AdaptiveSmoothRM(gt,  varargin)
    % smoothRateMap = SmoothRM(gt,  varargin)
    % adaptive smoothing - Skagg's 1996, Hippocampus
    if isempty(gt.pfObject), gt.LoadPF; end
    defClus = gt.pfObject.acceptedUnits;
    [type, IF_REPORTFIG, cluIdx, alpha, nBins, state, smoothFactor] = ...
        DefaultArgs(varargin, {'load', 0, defClus, 1e1, 50, 'RUN', 0.03});
    switch type
      case 'load'
        if FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename,'.mat'])
            load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename,'.mat']);
        end
        smoothRateMap = smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, cluIdx));
      case 'compute'

        if isempty(gt.pfObject)
            gt.LoadPF;
        end
        if isempty(gt.res)
            gt.LoadCR;
        end
        if isempty(cluIdx)
            cluIdx = gt.pfObject.acceptedUnits;
        else
            cluIdx = gt.pfObject.acceptedUnits(ismember(gt.pfObject.acceptedUnits, cluIdx));
            if isempty(cluIdx), return; end;
        end

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
            statePeriods = round(statePeriods .* gt.sampleRate ./ gt.lfpSampleRate) + 1;
            trialStartTime_pos = 0;
            markerNo = 1;
          case  'kenji'
            statePeriods = load([gt.paths.data, gt.filebase '.sts.', state]); % @lfp fs
            statePeriods = IntersectRanges(statePeriods, gt.trialPeriods);
            posStatePeriods = round(statePeriods .* gt.trackingSampleRate  ./ gt.lfpSampleRate) + 1;
            statePeriods = round(statePeriods .* gt.sampleRate ./ gt.lfpSampleRate) + 1;       
            trialStartTime_pos =  round(gt.trialPeriods(1, 1) .*  gt.trackingSampleRate ./ gt.lfpSampleRate) + 1;
            posStatePeriods = IntersectRanges(posStatePeriods, gt.goodPosPeriods + trialStartTime_pos); 
            markerNo = 1;
        end
        pos = SelectPeriods(sq(gt.position(:,markerNo,:)), posStatePeriods - trialStartTime_pos, 'c');
        %convert spike times to to tracking sample rate
        res = round(gt.res(ismember(gt.clu, cluIdx)) .* gt.trackingSampleRate ./ gt.sampleRate) + 1;
        clu = gt.clu(ismember(gt.clu, cluIdx));
        [res, resIdx] = SelectPeriods(res, posStatePeriods, 'd', 1, 1);
        clu = clu(resIdx);
        nClu = length(cluIdx);
        for lClu = 1 : nClu
            [~, d1, d2, d3, spkCnt(:, :, lClu), occCount(:, :, lClu)] = GenericPF.ComputeRateMap(gt, res(clu == cluIdx(lClu)), pos, [], .03, [], 0); 
        end
        matlabpool local 4
        parfor kClu = 1 : nClu
            if ~isempty(gt.pfObject.rateMap{cluIdx(kClu)})
                smoothMap = nan(length(gt.pfObject.xBin), length(gt.pfObject.yBin));
                for x = 1 : length(gt.pfObject.xBin)
                    for y = 1 : length(gt.pfObject.yBin)
                        radius = 1;
                        occSamples = GetPxInRadius(occCount(:, :, kClu) , x, y, radius);
                        nOccSamples = sum(~isnan(occSamples));
                        xySpkCnt =  sum(GetPxInRadius(spkCnt(:,:,kClu), x, y, radius));
                        condition = (alpha / (nOccSamples * sqrt(xySpkCnt)));
                        % expand the radius until r > condition
                        while radius <= condition | isinf(condition) | isnan(condition)
                            radius = radius + 1;
                            if radius > length(gt.pfObject.xBin) | radius > length(gt.pfObject.yBin), break; end
                            occSamples = GetPxInRadius(occCount, x, y, radius);
                            nOccSamples = sum(~isnan(occSamples));
                            xySpkCnt =  sum(GetPxInRadius(spkCnt(:,:,kClu), x, y, radius));
                            condition = (alpha / (nOccSamples * sqrt(xySpkCnt)));
                        end
                        smoothMap(x, y) = gt.trackingSampleRate * xySpkCnt / (nOccSamples * radius);                
                    end
                end                   
                smoothRateMap(:, :, kClu) = SmoothSurface(smoothMap, smoothFactor);
            end
        end
        matlabpool close
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename,'.mat'], 'smoothRateMap');
    end
end

