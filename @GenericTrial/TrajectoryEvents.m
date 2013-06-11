function trajEvntPeriods = TrajectoryEvents(gt, varargin)
% out = TrajectoryEvents(gt, varargin)
    if isempty(gt.pfObject), gt = gt.LoadPF; end
    defClus = gt.pfObject.acceptedUnits(gt.pfObject.sparsity < .35);
    [state, clus, binSize, bandWidth] = DefaultArgs(varargin, {'SWS', defClus, 50e-3, 30e-3});

    [res, clu] = gt.LoadStateRes(state);
    res = res(ismember(clu, clus));
    clu = clu(ismember(clu, clus));
    [spikeDensity, binEdges] = SpikeDensity(res, gt.sampleRate, binSize, bandWidth);
    mu = mean(spikeDensity);
    sigma = std(spikeDensity);
    trajEvnt = InOut(spikeDensity >= mu + 3 * sigma); % detect traj events
    trajEvnt = trajEvnt{1};
    aboveMeanPeriods = InOut(spikeDensity > mu);
    aboveMeanPeriods = aboveMeanPeriods{1};
    trajEvntBins = aboveMeanPeriods(RangeInInterval(trajEvnt, aboveMeanPeriods), :);
    trajEvntPeriods = round([binEdges(trajEvntBins(:, 1))', binEdges(trajEvntBins(:, 2))'] * gt.sampleRate);
end



