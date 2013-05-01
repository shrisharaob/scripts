function smoothRateMap = SmoothRM(gt,  varargin)
% Skagg's 1996, Hippocampus

    [cluIdx, alpha, nBins, state] = DefaultArgs(varargin, {[], 1e4, 50, 'RUN'});
    if isempty(gt.pfObject)
        gt = gt.LoadPF;
    end
    if isempty(gt.res)
        gt = gt.LoadCR;
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
        keyboard;
        statePeriods = load([gt.paths.data, gt.filebase '.sts.', state]); % @lfp fs
        statePeriods = IntersectRanges(statePeriods, gt.trialPeriods);
        posStatePeriods = round(statePeriods .* gt.trackingSampleRate  ./ gt.lfpSampleRate) + 1;
        statePeriods = round(statePeriods .* gt.sampleRate ./ gt.lfpSampleRate) + 1;       
        %posStatePeriods = round(statePeriods .* gt.trackingSampleRate  ./ gt.lfpSampleRate) + 1;
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

    %    [res, resIdx] = SelectPeriods(gt.res, gt.trialPeriods, 'd', 1, 1);
    % clu = gt.clu(resIdx);
    % [res, resIdx] = SelectPeriods(res, gt.goodPosPeriods, 'd', 1, 1);
    % clu = clu(resIdx);
    % pos = SelectPeriods(gt.position( :, 1 ,:), gt.goodPosPeriods, 'c');
    %    [~, occCount] = Occupancy(gt);
    for lClu = 1 : nClu
         [~, d1, d2, d3, spkCnt(:, :, lClu), occCount(:, :, lClu)] = GenericPF.ComputeRateMap(gt, res(clu == cluIdx(lClu)), pos, [], .03, [], 0); 
    end
    for kClu = 1 : nClu
        if ~isempty(gt.pfObject.rateMap{cluIdx(kClu)})
         
            for x = 1 : length(gt.pfObject.xBin)
                for y = 1 : length(gt.pfObject.yBin)
                    radius = 1;
                    occSamples = GetPxInRadius(occCount(:, :, kClu) , x, y, radius);
                    nOccSamples = sum(~isnan(occSamples));
                    xySpkCnt =  sum(GetPxInRadius(spkCnt(:,:,kClu), x, y, radius));
                    while radius > alpha / (nOccSamples * sqrt(xySpkCnt) + eps) 
                        radius = radius + 1;
                        occSamples = GetPxInRadius(occCount, x, y, radius);
                        nOccSamples = sum(~isnan(occSamples));
                        xySpkCnt =  sum(GetPxInRadius(spkCnt(:,:,kClu), x, y, radius));
                    end
                    smoothMap(x, y) = gt.trackingSampleRate * xySpkCnt / nOccSamples;                
                end
            end                   
            keyboard;
            smoothedRateMap(:, :, kClu) = smoothMap;
        end
    end
end

