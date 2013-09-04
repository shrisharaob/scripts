function [trajEvntPeriods, pars] = TrajectoryEvents(gt, IF_COMPUTE, preOrPost, varargin)
    % trajEvntPeriods = TrajectoryEvents(gt, IF_COMPUTE, preOrPost, varargin) @sample fs
    % returns periods with population events defined as a minimum
    % number of cells firing within a time window; preceeded and
    % followed by a minimum interval of silence
    % -------
    % Inputs:
    %    IF_COMPUTE
    %    preOrPost
    %    state
    %    clus2Select
    %    minEvntCells
    %    binSize
    %    overlap
    %     minSilenceWin  

    fprintf('\n traj evnts ...')
    trajEvntPeriods = [];
    pars = [];
    if isempty(gt.pfObject), gt.LoadPF; end
    if IF_COMPUTE, if isempty(gt.clu), gt.LoadCR; end; end
    defClus2Select = gt.pfObject.acceptedUnits(gt.pfObject.sparsity < 0.7);
    [state, clus2Select, minEvntCells, binSize, overlap, minSilenceWin] = ...
        DefaultArgs(varargin, {'SWS', defClus2Select, 5, 150e-3, 0, 50e-3 });

    if ~FileExists([gt.paths.analysis, gt.filebase, '.', state, '.' gt.trialName, '.', preOrPost, '.', mfilename, '.mat']), IF_COMPUTE = 1; end
    if ~IF_COMPUTE
        load([gt.paths.analysis, gt.filebase, '.', state, '.' gt.trialName, '.', preOrPost, '.', mfilename, '.mat']);
        fprintf('done ! \n');
        return;
    end
    sts  = gt.SleepPeriods(gt.sampleRate);
    switch preOrPost
      case 'pre'
        sts = sts{1};
      case 'post'
        sts = sts{2};
    end
    [res, resIdx] = SelectPeriods(gt.res, sts, 'd');
    clu  = gt.clu(resIdx);
    [res, resIdx] = sort(res);
    clu = clu(resIdx);
    res = res(ismember(clu, clus2Select));
    clu = clu(ismember(clu, clus2Select));
    if length(res) > 10
    res = res ./ gt.sampleRate; % res in sec
    binEdges = [[res(1) : binSize * (1 - overlap) : res(end)]',  [res(1) + binSize : binSize * (1 - overlap) : res(end) + binSize]'];
    nBins = size(binEdges, 1) - 1;
    IS_EVNT_BIN = false(nBins, 1);
    matlabpool local 8
    parfor mBin = 2 : nBins - 1 % the 1st and last bins cannot be tested for silence
    % candidate evnts identified by number of cells active in the bin
        resInBinIdx = res >= binEdges(mBin, 1) & res < binEdges(mBin, 2);
        counts(mBin) = sum(resInBinIdx);
        IS_EVNT_BIN(mBin) = length(unique(clu(resInBinIdx))) > minEvntCells;
        if IS_EVNT_BIN(mBin) % check if it is flanked by silence in minSilenceWin
            [preCount, ~] = histc(res, [binEdges(mBin, 1) - minSilenceWin, binEdges(mBin, 1)]);
            [sucCount, ~] = histc(res, [binEdges(mBin, 2), binEdges(mBin, 2) + minSilenceWin]);
            if preCount(1) >  1 | sucCount(1) > 1, IS_EVNT_BIN(mBin) = false; end 
        end
    end
    matlabpool close;
    evntBins = find(IS_EVNT_BIN);
    counts = counts ./ (nBins * binSize);
    trajEvntPeriods = round(binEdges(evntBins, :) .* gt.sampleRate);
    pars.clus2Select = clus2Select;
    pars.minEvntCells  = minEvntCells;
    pars.binSize = binSize;
    save([gt.paths.analysis, gt.filebase, '.', state, '.' gt.trialName, '.', preOrPost, '.', mfilename, '.mat'], 'pars', 'trajEvntPeriods');
    end
end




