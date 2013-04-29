function smoothRateMap = SmoothRM(gt,  varargin)
% Saki's 1996, Hippocampus

    [cluId, alpha, nBins] = DefaultArgs(varargin, {[], [], 50});
    if isempty(gt.pfObject)
        gt = gt.LoadPF;
    end
    if isempty(gt.res)
        fprintf('\n loading clures ...');
        gt = gt.LoadCR;
        fprintf(' done');
    end
    cluIdx = gt.pfObject.acceptedUnits;
    nClu = length(cluIdx);
    [res, resIdx] = SelectPeriods(gt.res, gt.trialPeriods, 'd', 1, 1);
    clu = gt.clu(resIdx);
    [res, resIdx] = SelectPeriods(res, gt.goodPosPeriods, 'd', 1, 1);
    clu = clu(resIdx);
    [~, occCount] = Occupancy(gt);
    for kClu = 1 : nClu
        if ~isempty(gt.pfObject.rateMap{cluIdx(kClu)})
            radius = 1;
            for x = 1 : length(gpf.xBins)
                for y = 1 : length(gpf.yBins)
                    occSamples = GetPxInRadius(occCount, x, y, radius);
                    nOccSamples = sum(~isnan(occSamples));
                    spkCnt = SpkCntAtPos(gt, res, ones(1, length(res)), [gt.pfObject.xBins(x), gt.pfObject.yBins(y)]);
                    while radius < alpha / (nOccSamples * sqrt(spkCnt)) 
                        occSamples = GetPxInRadius(occCount, x, y, radius);
                        nOccSamples = sum(~isnan(occSamples));
                        spkCnt = SpkCntAtPos(gt, res, ones(1, length(res)), [gt.pfObject.xBins(x), gt.pfObject.yBins(y)]);
                        radius = radius + 1;
                    end
                    smoothMap(x, y) = gt.trackinSampleRate * spkCnt / nOccSamples;                
                end
            end                   
            smoothedRateMap(:, :, kClu) = 
        end
    end
end

