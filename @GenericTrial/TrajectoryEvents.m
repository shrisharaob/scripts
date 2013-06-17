function trajEvntPeriods = TrajectoryEvents(gt, IF_COMPUTE, varargin)
% trajEvntPeriods = TrajectoryEvents(gt, varargin) @sample fs
    if ~IF_COMPUTE
        load([gt.paths.anslysis, gt.filebase, '.', state, '.' gt.trialName, '.', mfilename, '.mat']);
    end
    if isempty(gt.pfObject), gt.LoadPF; end
    defClus2Select = gt.pfObject.acceptedUnits(gt.pfObject.sparsity < .35);
    [state, clus2Select, minEvntCells, binSize, bandWidth] = ...
        DefaultArgs(varargin, {'SWS', defClus2Select, 5, 250e-3, 30e-3});

    [res, clu] = gt.LoadStateRes(state, [],[],[],1);
    [res, resIdx] = sort(res);
    clu = clu(resIdx);
    tres = res;
    tclu = clu;
%     res = tres(1 : 1e4);
%     clu = tclu(1 : 1e4);
    res = res(ismember(clu, clus2Select));
    clu = clu(ismember(clu, clus2Select));
    %    [spikeDensity, binEdges] = SpikeDensity(res, gt.sampleRate, binSize, bandWidth);
    res = res ./ gt.sampleRate;
    binEdges = res(1) : binSize : res(end) + binSize;
    nBins = length(binEdges) - 1;
    [counts, cntIdx] = histc(res, binEdges);
    for mBin = 1 : nBins
        % candidate evnts identified by number of cells active in the bin
        IS_EVNT_BIN(mBin) = length(unique(clu(cntIdx == mBin))) > minEvntCells; 
    end
    cellOrder = SortCellLoc(gt, clus2Select);
    evntBins = find(IS_EVNT_BIN);
    counts = counts ./ (nBins * binSize);
    trajEvntPeriods = round([binEdges(evntBins)', binEdges(evntBins + 1)'] .* gt.sampleRate);
    
    pars.clus2Select = clus2Select;
    pars.minEvntCells  = minEvntCells;
    pars.binSize = binSize;
    save([gt.paths.analysis, gt.filebase, '.', state, '.' gt.trialName, '.', mfilename, '.mat'], 'pars', 'trajEvntPeriods');
    keyboard;
%     mu = mean(spikeDensity);
%     sigma = std(spikeDensity);
%     trajEvnt = InOut(spikeDensity >= mu + 3 * sigma); % detect traj events
%     trajEvnt = trajEvnt{1};
%     aboveMeanPeriods = InOut(spikeDensity > mu);
%     aboveMeanPeriods = aboveMeanPeriods{1};
%     trajEvntBins = aboveMeanPeriods(RangeInInterval(trajEvnt, aboveMeanPeriods), :);
%     trajEvntPeriods = round([binEdges(trajEvntBins(:, 1))', binEdges(trajEvntBins(:, 2))'] * gt.sampleRate);
end



