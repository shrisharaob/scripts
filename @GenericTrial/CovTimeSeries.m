function CovTimeSeries(trial, ThPh,  varargin)
% CovTimeSeries(trial, ThPh,  varargin)
% [IF_REPORT_FIG, commonClus, timeWinSiz, overlap, pairs2analyse, IF_SAVE, thBinSiz, ccgBinSiz, fileTag, nThCycles]
    
    if nargin<1, help CovTimeSeries; return; end
    if isempty(trial.pfObject)
        trial = trial.LoadPF;
    end
    [IF_REPORT_FIG, commonClus, timeWinSiz, overlap, pairs2analyse, IF_SAVE, thBinSiz, ccgBinSiz, fileTag, nThCycles] = ...
        DefaultArgs(varargin, {0, [], 2000e-3, 0.8,[1:size(trial.pfObject.selectedPairs, 1)], 1, 10, 10e-3, [], 10});
    pfObject = trial.pfObject;
%     if FileExists(['~/data/analysis/' trial.filebase '/' trial.filebase '.SelectedPyrCells.mat'])
%         load(['~/data/analysis/' trial.filebase '/' trial.filebase '.SelectedPyrCells.mat']);
%         cellPairs = pfObject.selectedPairs(pairs2analyse, :);
%     else 
%         cellPairs =pfObject.selectedPairs;
%     end
%     switch trial.datasetType
%         case 'kenji'
%           filetag = GenFiletag(roi, arena);
%           if isempty(commonClus)
%               load(['~/data/analysis/kenji/', filebase, '/', filebase, filetag 'commonClus.mat']);
%           end
          if  isempty(commonClus) || length(commonClus) == 1, return; end
      % case  'MTA'
%         % commonClus        
%     end
    cellPairs = nchoosek(commonClus, 2);
    nPairs = size(cellPairs, 1);
    if isempty(trial.res)
        trial = trial.LoadCR;  
    end
    res = trial.res(ismember(trial.clu, commonClus)); % load res only for the units in roi
    clu = trial.clu(ismember(trial.clu, commonClus));
    res = round(res .* trial.lfpSampleRate ./ trial.sampleRate) + 1; % convert res to lfp sample rate
    [res, resIdx] = SelectPeriods(res, trial.trialPeriods, 'd', 1, 1);
    clu = clu(resIdx);
    thetaPeriods = trial.ThetaBoundaries(ThPh, commonClus, [], [], nThCycles);
    nCycles = size(thetaPeriods, 1);
    cv = zeros(nPairs, nCycles);
    oldStr = [];
    for kCv = 1 : nCycles   
        str = sprintf(['#', num2str(kCv) ' of ' num2str(nCycles)]);
        fprintf([repmat('\b', 1, length(oldStr)), str]);
        oldStr = str;
        [curRes, curResId] = SelectPeriods(res(ismember(clu, commonClus)), thetaPeriods(kCv, :),'d',1,1);
        curClu = clu(curResId);
        if ~isempty(curRes) & length(curRes) > 50
            for lPair = 1 : nPairs
                pRes = [];
                pClu = [];
                for kClu = 1 : 2
                    pRes = [pRes; curRes(ismember(curClu, cellPairs(lPair, kClu)))];
                    pClu = [pClu; cellPairs(lPair, kClu) * ones(length( curRes(ismember(curClu, cellPairs(lPair, kClu)))) , 1) ];
                end
                if length(pRes) > 20
                    halfBins = round(ccgBinSiz * 1e1 / ccgBinSiz);
                    nBinSamples = round(ccgBinSiz * trial.sampleRate);
                    [ccgOut(:, :, :, kCv), ccgTimeAx, Pairs] = ...
                        myCCG(pRes, pClu, nBinSamples, halfBins, trial.sampleRate, cellPairs(lPair, :) , 'count');
                end
                keyboard;  
            end 
            keyboard;
        end
    end
end
%     timeWinSiz = round(timeWinSiz * trial.sampleRate);
%     startEdges = trial.trialPeriods(1) : timeWinSiz * (1 - overlap) : trial.trialPeriods(2) - timeWinSiz;
%     endEdges = trial.trialPeriods(1) + timeWinSiz : timeWinSiz * (1 - overlap) : trial.trialPeriods(2);
%     winEdges = [startEdges', endEdges'];
%     nWindows = size(winEdges, 1); % total no of bins
%     for kTimeWin = 1 : nWindows
%         [kWinRes, kWinResIdx] = SelectPeriods(res, winEdges(kTimeWin, :), 'd', 1, 1);
%         kWinClu = clu(kWinResIdx);
%         for lPair = 1 : nPairs
%             pRes = [];
%             pClu = [];
%             for kClu = 1 : 2
%                 pRes = [pRes; kWinRes(kWinClu == cellPairs(lPair, kClu))];
%                 pClu = [pClu; cellPairs(lPair, kClu) * ones(length( kWinRes(kWinClu == cellPairs(lPair, kClu))), 1)];
%             end
%             %            if ~isempty(pClu), keyboard; end
%             halfBins = round(ccgBinSiz * 1e1 / ccgBinSiz);
%             % if isinf(halfBins), keyboard; end
%             nBinSamples = round(ccgBinSiz * trial.sampleRate);
%             [ccgOut, ccgTimeAx, Pairs] = ...
%                 myCCG(pRes, pClu, nBinSamples, halfBins, trial.sampleRate, cellPairs(lPair, :) , 'count');
%         end 
%     end




