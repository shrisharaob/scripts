function CovTimeSeries(trial, varargin)
% CovTimeSeries(trial, varargin)
% [ThPh, commonClus, IF_REPORT_FIG, nThCycles, overlap, pairs2analyse, IF_SAVE, binSize, fileTag] 
% [], [], 0, 3, 0.8,[1:size(trial.pfObject.selectedPairs, 1)], 1, 10e-3, []
% computes cov in theta cycles  
     
    if nargin<1, help CCGofSegments; return; end
    if isempty(trial.pfObject), trial = trial.LoadPF; end
    [ThPh, commonClus, IF_REPORT_FIG, nThCycles, overlap, pairs2analyse, IF_SAVE, binSize, fileTag] = ...
        DefaultArgs(varargin, {[], [], 0, 3, 0.8,[1:size(trial.pfObject.selectedPairs, 1)], 1, 10e-3, []});

    pfObject = trial.pfObject;
    if FileExists(['~/data/analysis/' trial.filebase '/' trial.filebase '.SelectedPyrCells.mat'])
        load(['~/data/analysis/' trial.filebase '/' trial.filebase '.SelectedPyrCells.mat']);
        cellPairs = pfObject.selectedPairs(pairs2analyse, :);
    else 
        cellPairs =pfObject.selectedPairs;
    end
    
    nPairs = size(cellPairs, 1);
    if isempty(trial.res), trial = trial.LoadCR; end
    if isempty(ThPh), load([trial.paths.data, trial.filebase, '.thpar.mat']);
    [res, resIdx] = SelectPeriods(trial.res, round(trial.trialPeriods .* trial.sampleRate ./ trial.lfpSampleRate) + 1 , 'd', 1, 1);
    clu = trial.clu(resIdx);
    winEdges = ConvertFs(trial.ThetaBoundaries(ThPh, commonClus, [], [], nThCycles), trial.lfpSampleRate, trial.sampleRate);
%   timeWinSiz = round(timeWinSiz * trial.sampleRate);
%   startEdges = trial.trialPeriods(1) : timeWinSiz * (1 - overlap) : trial.trialPeriods(2) - timeWinSiz;
%   endEdges = trial.trialPeriods(1) + timeWinSiz : timeWinSiz * (1 - overlap) : trial.trialPeriods(2);
%   winEdges = [startEdges', endEdges'];
    nWindows = size(winEdges, 1); % total no of bins
    oldStr = [];
    count = 0;
    for kTimeWin = 1 : nWindows
        [kWinRes, kWinResIdx] = SelectPeriods(res, winEdges(kTimeWin, :), 'd');
        kWinClu = clu(kWinResIdx);
        str = sprintf(['#', num2str(kTimeWin) ' of ' num2str(nWindows)]);
        fprintf([repmat('\b', 1, length(oldStr)), str]);
        oldStr = str;
        [binnedPos, coverageMask] = Coverage(gt, roi, arena, markerNo, 0, 0);     
        for lPair = 1 : nPairs         
            pRes = [];
            pClu = [];
            for kClu = 1 : 2
                pRes = [pRes; kWinRes(kWinClu == cellPairs(lPair, kClu))];
                pClu = [pClu; cellPairs(lPair, kClu) * ones(length( kWinRes(kWinClu == cellPairs(lPair, kClu))), 1)];
            end
            if sum(pClu == cellPairs(lPair, 1)) > 10 && sum(pClu == cellPairs(lPair, 2)) > 10
                halfBins = round(binSize * 1e1 / binSize);
                binSize = round(binSize * trial.sampleRate);
                try
                    [ccgOut, ccgTimeAx, Pairs] =  myCCG(pRes, pClu, binSize, halfBins, trial.sampleRate, cellPairs(lPair, :), 'count');
                    ccg(:, kTimeWin) = ccgOut(:, 1, 2);
                    count = count + 1;
                catch err
                    
                end
            end
        end        
    end
    keyboard;
end
