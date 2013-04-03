    % load ([gt.paths.data, gt.filebase '.thpar.mat'] )
    thetaPhase = SelectPeriods(ThPh, gt.goodPosPeriods, 'c', 1); % theta phase in lfp sample rate
    binSize = 10; % 10 deg
    nBins = 360 / binSize;
    bins = linspace(0, pi, nBins);
    [~, binIdx] = histc(thetaPhase, bins);
    if isempty(gt.res)
        [gt.res, gt.clu] = LoadCluRes([gt.paths.data, gt.filebase]);
    end
    res = SelectPeriods(gt.res, gt.trialPeriods, 'd', [], 1);
    res = round(res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
    res = 