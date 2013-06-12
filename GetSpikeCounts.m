function [sc, winEdges, varargout] = GetSpikeCounts(trial, varargin)
%  returns spike counts in wins
% trial - generic trial class object
% winSize in secs
% statePeriod - in samplerate
% overlap - overlap factor of windows
% function sc = GetSpikeCounts(trial, winSize, cluIdx)
    format long;
    switch trial.datasetType
      case 'kenji'
        defMarker = 1;
      case 'MTA'
        defMarker = 7;
    end
    [winSiz, statePeriods, cluIdx, overlap, markerNo] =  ...
        DefaultArgs(varargin, {200e-3, [], [], 0, defMarker});
    %  xySegLength = size(trial.position, 1);
    %  winSiz = round(winSiz * trial.sampleRate);
    %  nWins = ceil(xySegLength ./ (trial.trackingSampleRate * winSize));
    %  winEdges = linspace(0, xySegLength / trial.trackingSampleRate, nWins);
    nClus = length(cluIdx);
    if isempty(trial.res)
        trial = trial.LoadCR;
    end
    [res, origIdx] = SelectPeriods(trial.res, statePeriods, 'd');
    xy = SelectPeriods(sq(trial.position(:, markerNo, :)), trial.goodPosPeriods, 'c');
    clu = trial.clu(origIdx);
    res = res(ismember(clu, cluIdx));
    clu = clu(ismember(clu, cluIdx));
    res = res ./ trial.sampleRate; % res in seconds
    startEdges = statePeriods(1, 1) ./ trial.sampleRate : winSiz * (1 - overlap) : statePeriods(end, 2) ./ trial.sampleRate - winSiz;
    endEdges = statePeriods(1, 1) ./ trial.sampleRate + winSiz : winSiz * (1 - overlap) : statePeriods(end, 2) ./ trial.sampleRate;
    winEdges = [startEdges', endEdges'];
    winEdgesATPosFs = round(RecenterPeriods(winEdges) .* trial.trackingSampleRate) + 1 ;
    nWins = size(winEdges, 1);
    sc = zeros(nClus, nWins);   
    for mWin = 1 : nWins
        for kClu = 1 : nClus
            kRes  = res(clu == cluIdx(kClu));
            % sc(kClu, mWin) = histc(kRes, winEdges(mWin, :));
            sc(kClu, mWin) = sum(kRes >= winEdges(mWin, 1)  & kRes <= winEdges(mWin, 2));
        end
        if winEdges(mWin, 2) <= size(xy, 1)
            posInWin(mWin, :) = nanmean(xy(winEdgesATPosFs(mWin, 1): winEdgesATPosFs(mWin, 2), :), 1);
        end
    end
    varargout{1} = posInWin;
end
