function out = PairReactivation(gt, varargin)
% out = PairReactivation(gt, varargin)
% returns pairwise reactivation likelihood

    out = [];
    [prePost, type, winSize, overlap ] = DefaultArgs(varargin, {'pre', 'load', 100e-3, 0});

    if isempty(gt.pfObject), gt.LoadPF; end
    switch type
      case 'compute'
        [sleepRes, sleepClu] = gt.LoadStateRes('SWS', [], [], [], 1, prePost);
        [runRes, runClu] = gt.LoadStateRes('RUN', 1, [], [], 1);
        
        cmnClus = intersect(unique(runClu), unique(sleepClu));
        if cmnClus > 1
            cellPairs = nchoosek(cmnClus, 2);
            % cellPairs = cellPairs(ismember(cellPairs, gt.pfObject.selectPairs));
            nPairs = size(cellPairs, 1);
        else
            fprintf('\n no cell pairs firing in both sleep and run \n');
            return;
        end

        %% common cells
        sleepRes = sleepRes(ismember(sleepClu, cmnClus));
        sleepClu = sleepClu(ismember(sleepClu, cmnClus));
        [sleepRes, resIdx] = sort(sleepRes);
        sleepClu = sleepClu(resIdx);
        
        runRes = runRes(ismember(runClu, cmnClus));
        runClu = runClu(ismember(runClu, cmnClus));
        [runRes, resIdx] = sort(runRes);
        runClu = runClu(resIdx);
        
        %%  population events
        %    [evntPeriods, ~] = gt.TrajectoryEvents(0, prePost, 'SWS');

        %% cofiring events in RUN
        runRes = runRes ./ gt.sampleRate; % res in sec
        binEdges = [[runRes(1) : winSize * (1 - overlap) : runRes(end)]',  [runRes(1) + winSize : winSize * (1 - overlap) : runRes(end) + winSize]'];
        nRunBins = size(binEdges, 1) - 1;
        IS_EVNT_BIN_RUN = false(nPairs, nRunBins);
        matlabpool local 8
        parfor kBin = 1 : nRunBins
            resInBinIdx = runRes >= binEdges(kBin, 1) & runRes < binEdges(kBin, 2);
            clusInBin = runClu(resInBinIdx);
            if clusInBin > 1
                pairsInBin = nchoosek(clusInBin, 2);
                IS_EVNT_BIN_RUN(:, kBin) = ismember(cellPairs, pairsInBin, 'rows');
            end
        end        
        matlabpool close
        
        %% cofiring events in sleep
        sleepRes = sleepRes ./ gt.sampleRate;
        binEdges = [[sleepRes(1) : winSize * (1 - overlap) : sleepRes(end)]',  [sleepRes(1) + winSize : winSize * (1 - overlap) : sleepRes(end) + winSize]'];
        nSleepBins = size(binEdges, 1) - 1;
        IS_EVNT_BIN_SLEEP = false(nPairs, nSleepBins);
        matlapool local 8
        parfor kBin = 1 : nSleepBins
            resInBinIdx = sleepRes >= binEdges(kBin, 1) & sleepRes < binEdges(kBin, 2);
            clusInBin = sleepClu(resInBinIdx);
            if clusInBin > 1
                pairsInBin = nchoosek(clusInBin, 2);
                IS_EVNT_BIN_SLEEP(:, kBin) = ismember(cellPairs, pairsInBin, 'rows');
            end
        end        
        matlabpool close

        %% likelihood of cofiring
        cofiringProb = (sum(IS_EVNT_BIN_RUN, 2) + sum(IS_EVNT_BIN_SLEEP, 2))  ./ (nRunBins + nSleepBins); 
        out.cellPairs = cellPairs;
        out.cofiringProb = cofiringProb;
        %out.evntBins = 
        keyboard;
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat'], 'out');

      case 'load'
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat']);
        
    end
end