function out = PairReactivation(gt, varargin)
% out = PairReactivation(gt, varargin)
% returns pairwise reactivation likelihood

    out = [];
    [prePost, type, nResample, IF_PLOT, winSize, overlap ] = ...
        DefaultArgs(varargin, {'pre', 'load', 0, false, 100e-3, 0});

    if isempty(gt.pfObject), gt.LoadPF; end
    if ~FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat']), type = 'compute'; end
    switch type
      case 'compute'
        [sleepRes, sleepClu] = gt.LoadStateRes('SWS', [], [], [], 1, prePost);
        runClu = gt.pfObject.acceptedUnits(gt.pfObject.sparsity < 0.4);
        cmnClus = intersect(unique(runClu), unique(sleepClu));
        if length(cmnClus) > 1
            cellPairs = nchoosek(cmnClus, 2);
            nPairs = size(cellPairs, 1);
        else
            fprintf('\n no cell pairs firing in both sleep and run \n');
            save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat'], 'out');
            return;
        end
        %% common cells
        sleepRes = sleepRes(ismember(sleepClu, cmnClus));
        sleepClu = sleepClu(ismember(sleepClu, cmnClus));
        [sleepRes, resIdx] = sort(sleepRes);
        sleepClu = sleepClu(resIdx);
        %% cofiring events in sleep
        sleepRes = sleepRes ./ gt.sampleRate;
        %          %% TEST
        %         sleepRes = sleepRes(1:1e3);
        %         sleepClu = sleepClu(1:1e3);
        
        binEdges = [[sleepRes(1) : winSize * (1 - overlap) : sleepRes(end)]',  [sleepRes(1) + winSize : winSize * (1 - overlap) : sleepRes(end) + winSize]'];
        nSleepBins = size(binEdges, 1) - 1;
        IS_EVNT_BIN_SLEEP = false(nPairs, nSleepBins);
        jitter = @(jitterWinSiz, x) x + (- jitterWinSiz  / 2) + jitterWinSiz .* rand(size(x)); 
        SPIKE_IN_BIN = false(length(cmnClus), nSleepBins);
        for kCell = 1 : length(cmnClus)
            kCellRes = sleepRes(sleepClu == cmnClus(kCell));
            for kBin = 1 : nSleepBins
                SPIKE_IN_BIN(kCell, kBin) = logical(sum(kCellRes >= binEdges(kBin, 1) & kCellRes < binEdges(kBin, 2)));
            end
        end
        for kPr = 1 : nPairs
            IS_EVNT_BIN_SLEEP(kPr, :) = SPIKE_IN_BIN(cmnClus == cellPairs(kPr, 1), :) & SPIKE_IN_BIN(cmnClus == cellPairs(kPr, 2), :);
        end
        if nResample
            chanceLvlProb = zeros(1,nPairs);
            fprintf('\n resampling ... \n');
            matlabpool local 8
            chanceLevelCnt = zeros(1, nSleepBins);
            parfor kResample = 1 : nResample
                SHUFFLED_SPIKE_IN_BIN = false(length(cmnClus), nSleepBins);
                SHUFFLED_EVNT_BIN = false(nPairs, nSleepBins);
                for kCell = 1 : length(cmnClus)
                    kCellRes = jitter(winSize, sleepRes(sleepClu == cmnClus(kCell)));
                    for kBin = 1 : nSleepBins
                        SHUFFLED_SPIKE_IN_BIN(kCell, kBin) = logical(sum(kCellRes >= binEdges(kBin, 1) & kCellRes < binEdges(kBin, 2)));
                    end
                end
                for kPr = 1 : nPairs
                    SHUFFLED_EVNT_BIN(kPr, :) =  SHUFFLED_SPIKE_IN_BIN(cmnClus == cellPairs(kPr, 1), :) & SHUFFLED_SPIKE_IN_BIN(cmnClus == cellPairs(kPr, 2), :);
                end
                chanceCnt(:, :, kResample) = SHUFFLED_EVNT_BIN;
            end
        end
        matlabpool close
        chanceLvlProb = sum(sum(chanceCnt, 3), 2)  ./ (nResample * nSleepBins);
        %% likelihood of cofiring
        out.cellPairs = cellPairs;
        out.chanceLvlProb = chanceLvlProb;
        out.dataCofiring = sum(IS_EVNT_BIN_SLEEP, 2) ./ length(IS_EVNT_BIN_SLEEP);
        out.cofiringAbvChance =  out.dataCofiring - out.chanceLvlProb;
        out.cofiringAbvChance(out.cofiringAbvChance < 0) = 0;
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat'], 'out');
      case 'load'
        if FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat'])
            load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat']);
        else, out = [];
        end
    end
end