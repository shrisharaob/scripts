function out = PairReactivation(gt, varargin)
% out = PairReactivation(gt, varargin)
% returns pairwise reactivation likelihood

    out = [];
    [prePost, type, IF_PLOT, winSize, overlap ] = DefaultArgs(varargin, {'pre', 'load', false, 100e-3, 0});

    if isempty(gt.pfObject), gt.LoadPF; end
    if ~FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat']), type = 'compute'; end
    switch type
      case 'compute'
        [sleepRes, sleepClu] = gt.LoadStateRes('SWS', [], [], [], 1, prePost);
        [runRes, runClu] = gt.LoadStateRes('RUN', 1, [], [], 1);
        
        cmnClus = intersect(unique(runClu), unique(sleepClu));
        if length(cmnClus) > 1
            cellPairs = nchoosek(cmnClus, 2);
            % cellPairs = cellPairs(ismember(cellPairs, gt.pfObject.selectPairs));
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
            if length(clusInBin) > 1
                pairsInBin = nchoosek(clusInBin, 2);
                %if size(pairsInBin, 2) == 1, keyboard; end
                IS_EVNT_BIN_RUN(:, kBin) = ismember(cellPairs, pairsInBin, 'rows');
            end
        end        
        matlabpool close
        
        %% cofiring events in sleep
        sleepRes = sleepRes ./ gt.sampleRate;
        binEdges = [[sleepRes(1) : winSize * (1 - overlap) : sleepRes(end)]',  [sleepRes(1) + winSize : winSize * (1 - overlap) : sleepRes(end) + winSize]'];
        nSleepBins = size(binEdges, 1) - 1;
        IS_EVNT_BIN_SLEEP = false(nPairs, nSleepBins);
        matlabpool local 8
        parfor kBin = 1 : nSleepBins
            resInBinIdx = sleepRes >= binEdges(kBin, 1) & sleepRes < binEdges(kBin, 2);
            clusInBin = sleepClu(resInBinIdx);
            if length(clusInBin) > 1
                pairsInBin = nchoosek(clusInBin, 2);
                IS_EVNT_BIN_SLEEP(:, kBin) = ismember(cellPairs, pairsInBin, 'rows');
            end
        end        
        matlabpool close
        %% likelihood of cofiring
        cofiringProb = (sum(IS_EVNT_BIN_RUN, 2) + sum(IS_EVNT_BIN_SLEEP, 2))  ./ (nRunBins + nSleepBins); 
        out.cellPairs = cellPairs;
        out.cofiringProb = cofiringProb;
        out.IS_EVNT_BIN_RUN = IS_EVNT_BIN_RUN;
        out.IS_EVNT_BIN_SLEEP = IS_EVNT_BIN_SLEEP;
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat'], 'out');
      case 'load'

        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.' mfilename, '.mat']);
        out.runCfProb = sum(out.IS_EVNT_BIN_RUN, 2) ./ length(out.IS_EVNT_BIN_RUN);
        out.sleepCfProb = sum(out.IS_EVNT_BIN_SLEEP, 2) ./ length(out.IS_EVNT_BIN_SLEEP);
        %        out.cofiringProb = 
        if IF_PLOT
            plot(out.runCfProb, out.sleepCfProb, 'ob', 'MarkerSize', 5);
            xlim([0, 1]);
            ylim([0, 1]);
            hold on;
            line;
            xlabel('RUN cofiring probablity');
            ylabel('SWS cofiring probablity');
            reportfig(gcf, mfilename, 0, [gt.filebase,  '    ', gt.trialName]);
            clf;
        end
    end
end