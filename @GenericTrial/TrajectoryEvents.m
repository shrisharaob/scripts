function trajEvntPeriods = TrajectoryEvents(gt, IF_COMPUTE, varargin)
% trajEvntPeriods = TrajectoryEvents(gt, varargin) @sample fs
    fprintf('\n traj evnts \n')
    if ~IF_COMPUTE
        load([gt.paths.anslysis, gt.filebase, '.', state, '.' gt.trialName, '.', mfilename, '.mat']);
    end
    if isempty(gt.pfObject), gt.LoadPF; end
    defClus2Select = gt.pfObject.acceptedUnits(gt.pfObject.sparsity < .5);
    [state, clus2Select, minEvntCells, binSize, minSilenceWin] = ...
        DefaultArgs(varargin, {'SWS', defClus2Select, 5, 250e-3, 50e-3});

    [res, clu] = gt.LoadStateRes(state, [],[],[],1);
    [res, resIdx] = sort(res);
    clu = clu(resIdx);
    tres = res(1 : 1e4);
    tclu = clu(1 : 1e4);
    res = res(ismember(clu, clus2Select));
    clu = clu(ismember(clu, clus2Select));
    %    [spikeDensity, binEdges] = SpikeDensity(res, gt.sampleRate, binSize, bandWidth);
    res = res ./ gt.sampleRate; % res in sec
    binEdges = res(1) : binSize : res(end) + binSize;
    nBins = length(binEdges) - 1;
    [counts, cntIdx] = histc(res, binEdges);
    for mBin = 2 : nBins - 1 % the 1st and last bins cannot be tested for silence
        % candidate evnts identified by number of cells active in the bin
        IS_EVNT_BIN(mBin) = length(unique(clu(cntIdx == mBin))) > minEvntCells; 
        if IS_EVNT_BIN(mBin) % check if it is flanked by silence in minSilenceWin
            [preCount, ~] = histc(res, [binEdges(mBin) - minSilenceWin, binEdges(mBin)]);
            [sucCount, ~] = histc(res, [binEdges(mBin), binEdges(mBin) + minSilenceWin]);
            if preCount > 1 & sucCount > 1, IS_EVNT_BIN(mBin) = false; end 
        end
    end


    cellOrder = SortCellLoc(gt, clus2Select);
    evntBins = find(IS_EVNT_BIN);
    counts = counts ./ (nBins * binSize);
    trajEvntPeriods = round([binEdges(evntBins)', binEdges(evntBins + 1)'] .* gt.sampleRate);
    
    pars.clus2Select = clus2Select;
    pars.minEvntCells  = minEvntCells;
    pars.binSize = binSize;
    save([gt.paths.analysis, gt.filebase, '.', state, '.' gt.trialName, '.', mfilename, '.mat'], 'pars', 'trajEvntPeriods');
    %     mu = mean(spikeDensity);
%     sigma = std(spikeDensity);
%     trajEvnt = InOut(spikeDensity >= mu + 3 * sigma); % detect traj events
%     trajEvnt = trajEvnt{1};
%     aboveMeanPeriods = InOut(spikeDensity > mu);
%     aboveMeanPeriods = aboveMeanPeriods{1};
%     trajEvntBins = aboveMeanPeriods(RangeInInterval(trajEvnt, aboveMeanPeriods), :);
%     trajEvntPeriods = round([binEdges(trajEvntBins(:, 1))', binEdges(trajEvntBins(:, 2))'] * gt.sampleRate);
end



