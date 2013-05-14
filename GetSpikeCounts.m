function [sc, winEdges] = GetSpikeCounts(trial, varargin)
%  returns spike counts in wins
% trial - generic trial class object
% winSize in secs
% statePeriod - in secs
% overlap - overlap factor of windows
% function sc = GetSpikeCounts(trial, winSize, cluIdx)
    
    [winSiz, statePeriods, cluIdx, overlap] =  ...
        DefaultArgs(varargin, {200e-3, [], [], 0});
    %  xySegLength = size(trial.position, 1);
    %  winSiz = round(winSiz * trial.sampleRate);
    %  nWins = ceil(xySegLength ./ (trial.trackingSampleRate * winSize));
    %  winEdges = linspace(0, xySegLength / trial.trackingSampleRate, nWins);
    nClus = length(cluIdx);
    if isempty(trial.res)
        trial = trial.LoadCR;
    end
    [res, origIdx] = SelectPeriods(trial.res, statePeriods, 'd');
    clu = trial.clu(origIdx);
    res = res(ismember(clu, cluIdx));
    clu = clu(ismember(clu, cluIdx));
    res = res ./ trial.sampleRate; % res in seconds
    startEdges = statePeriods(1, 1) ./ trial.sampleRate : winSiz * (1 - overlap) : statePeriods(end, 2) ./ trial.sampleRate - winSiz;
    endEdges = statePeriods(1, 1) ./ trial.sampleRate + winSiz : winSiz * (1 - overlap) : statePeriods(end, 2) ./ trial.sampleRate;
    winEdges = [startEdges', endEdges'];
    nWins = size(winEdges, 1);
    sc = zeros(nClus, nWins);   
    for mWin = 1 : nWins
        for kClu = 1 : nClus
            kRes  = res(clu == cluIdx(kClu));
            %            sc(kClu, mWin) = histc(kRes, winEdges(mWin, :));
            sc(kClu, mWin) = sum(kRes >= winEdges(mWin, 1)  & kRes <= winEdges(mWin, 2));
        end
    end
end
