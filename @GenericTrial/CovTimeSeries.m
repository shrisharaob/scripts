function CovTimeSeries(trial, varargin)
% 
     
    if nargin<1, help CCGofSegments; return; end
    if isempty(trial.pfObject)
        trial = trial.LoadPF;
    end
    [IF_REPORT_FIG, timeWinSiz, overlap, pairs2analyse, IF_SAVE, binSize, fileTag] = ...
        DefaultArgs(varargin, {0, 2000e-3, 0.8,[1:size(trial.pfObject.selectedPairs, 1)], 1, 10e-3, []});

    pfObject = trial.pfObject;
    if FileExists(['~/data/analysis/' trial.filebase '/' trial.filebase '.SelectedPyrCells.mat'])
        load(['~/data/analysis/' trial.filebase '/' trial.filebase '.SelectedPyrCells.mat']);
        cellPairs = pfObject.selectedPairs(pairs2analyse, :);
    else 
        cellPairs =pfObject.selectedPairs;
    end
    nPairs = size(cellPairs, 1);
    if isempty(trial.res)
        trial = trial.LoadCR;  
    end
    [res, resIdx] = SelectPeriods(trial.res, round(trial.trialPeriods .* trial.sampleRate ./ trial.lfpSampleRate) + 1 , 'd', 1, 1);
    clu = trial.clu(resIdx);
    timeWinSiz = round(timeWinSiz * trial.sampleRate);
    startEdges = trial.trialPeriods(1) : timeWinSiz * (1 - overlap) : trial.trialPeriods(2) - timeWinSiz;
    endEdges = trial.trialPeriods(1) + timeWinSiz : timeWinSiz * (1 - overlap) : trial.trialPeriods(2);
    winEdges = [startEdges', endEdges'];
    nWindows = size(winEdges, 1); % total no of bins
    for kTimeWin = 1 : nWindows
        [kWinRes, kWinResIdx] = SelectPeriods(res, winEdges(kTimeWin, :), 'd', 1, 1);
        kWinClu = clu(kWinResIdx);
        for lPair = 1 : nPairs
            pRes = [];
            pClu = [];
            for kClu = 1 : 2
                pRes = [pRes; kWinRes(kWinClu == cellPairs(lPair, kClu))];
                pClu = [pClu; cellPairs(lPair, kClu) * ones(length( kWinRes(kWinClu == cellPairs(lPair, kClu))), 1)];
            end
            %            if ~isempty(pClu), keyboard; end
            halfBins = round(binSize * 1e1 / binSize);
            % if isinf(halfBins), keyboard; end
            nBinSamples = round(binSize * trial.sampleRate);
            [ccgOut, ccgTimeAx, Pairs] = ...
                myCCG(pRes, pClu, nBinSamples, halfBins, trial.sampleRate, cellPairs(lPair, :) , 'count');
        end 
        keyboard;
    end
keyboard;
end
