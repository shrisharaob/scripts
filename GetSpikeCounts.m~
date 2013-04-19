function [sc, binCenters] = GetSpikeCounts(trial, varargin)
    %  returns spike counts in bins
    % trial - generic trial class object
    % binSize in secs
    % function sc = GetSpikeCounts(trial, binSize, cluIdx)
    
    [binSize, cluIdx, state] = DefaultArgs(varargin, {500e-3, [], 'walk'});
    xySegLength = size(trial.position, 1);
    nBins = ceil(xySegLength ./ (trial.trackingSampleRate * binSize));
    binCenters = linspace(0, xySegLength / trial.trackingSampleRate, nBins);
    nClus = length(cluIdx);
    sc = zeros(nClus, nBins);
    for i = 1:length(trial.states),
        if strcmp(trial.states{i}.label,state),
            statePeriods = trial.states{i}.statePeriods;
            break;
        end
    end
    %keyboard;
    statePeriods = statePeriods + 1; % bug -  xyz index starts from 0
    if strcmp(trial.datasetType, 'MTA')
        statePeriods = round(statePeriods .* (trial.sampleRate / ...
        trial.trackingSampleRate))+ 1;
    else % state periods in lfp sampleing rate
        statePeriods = round(statePeriods .* (trial.sampleRate / ...
        trial.lfpSampleRate))+ 1;
    end
    
    [res, origIdx] = SelectPeriods(trial.res, statePeriods);
    res = res ./ trial.sampleRate; % res in seconds
    clu = trial.clu(origIdx);
    for kClu = 1 : nClus
        kRes  = res(clu == cluIdx(kClu));
        sc(kClu, :) = hist(kRes,binCenters);
    end
end





