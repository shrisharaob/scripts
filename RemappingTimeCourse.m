
    load ([gt.paths.data, gt.filebase '.thpar.mat'] );

    
    if isempty(gt.res)
        [gt.res, gt.clu] = LoadCluRes([gt.paths.data, gt.filebase]);
    end
    %   thetaPhase = SelectPeriods(ThPh, gt.trialPeriods, 'c'); % theta phase in lfp sample rate
    
    res = round(gt.res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
res = SelectPeriods(gt.res, gt.trialPeriods, 'd');
res = SelectPeriods(res, round(gt.goodPosPeriods .* gt.lfpSampleRate ...
                               ./ gt.trackingSampleRate) + 1 + ...
                    res(1),'d');
stateThPh = ThPh(res);

    binSize = 10; % 10 deg
    nBins = 360 / binSize;
    bins = linspace(0, pi, nBins);
    [count, binIdx] = histc(stateThPh, bins);
		
