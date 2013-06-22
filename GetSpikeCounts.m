function [sc, winEdges, varargout] = GetSpikeCounts(trial, varargin)
% function sc = GetSpikeCounts(trial, varargin) 
% [winSiz, statePeriods, cluIdx, overlap, markerNo] 
% {200e-3, [], [], 0, defMarker}
%  returns spike counts in wins
% trial - generic trial class object
% winSize in secs
% statePeriod - in samplerate
% overlap - overlap factor of windows

    format long;
    switch trial.datasetType
      case 'kenji'
        defMarker = 1;
      case 'MTA'
        defMarker = 7;
    end
    [winSiz, state, cluIdx, overlap, markerNo] =  ...
        DefaultArgs(varargin, {200e-3, [], [], 0, defMarker});
    nClus = length(cluIdx);
    if isempty(trial.res)
        trial = trial.LoadCR;
    end
    [res, clu, xy] = trial.LoadStateRes(state, 1, 0);
    res = res(ismember(clu, cluIdx));
    clu = clu(ismember(clu, cluIdx));
    res = res ./ trial.sampleRate; % res in seconds
    edges = 1 : winSiz * (1 - overlap) : res(end) + winSiz;
    winEdges = [edges(1: end-1)', edges(2 :end)'];
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
        else 
            posInWin(mWin, :) = [nan, nan];
        end
    end
    varargout{1} = posInWin;
end
