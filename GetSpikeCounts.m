function [sc, binEdges] = GetSpikeCounts(trial, varargin)
    %  returns spike counts in bins
    % trial - generic trial class object
    % binSize in secs
    % function sc = GetSpikeCounts(trial, binSize, cluIdx)
    
    [binSize, cluIdx, state] = DefaultArgs(varargin, {200e-3, [], 'RUN'});
    xySegLength = size(trial.position, 1);
    nBins = ceil(xySegLength ./ (trial.trackingSampleRate * binSize));
    binEdges = linspace(0, xySegLength / trial.trackingSampleRate, nBins);
    nClus = length(cluIdx);
    sc = zeros(nClus, nBins);
    switch trial.datasetType
      case 'MTA'
        for i = 1:length(trial.states)
            if strcmp(state, 'RUN')
                state = 'walk';
            end
            if strcmp(trial.states{i}.label, state),
                statePeriods = trial.states{i}.statePeriods;
                break;
            end
        end
      case 'kenji'
        for i = 1:length(trial.states)   
            if strcmp(trial.states{i}.name, state)
                statePeriods = trial.states{i}.statePeriods; % @ lfp fs
                break;
            end
        end
    end
    if strcmp(trial.datasetType, 'MTA')
        statePeriods = round(statePeriods .* (trial.sampleRate /  trial.trackingSampleRate))+ 1;
    else % state periods in lfp sampling rate
        statePeriods = round(statePeriods .* (trial.sampleRate /  trial.lfpSampleRate))+ 1;  
    end
    if isempty(trial.res)
       trial = trial.LoadCR;
    end
    [res, origIdx] = SelectPeriods(trial.res, statePeriods, 'd', 1, 1);
    res = res ./ trial.sampleRate; % res in seconds
    clu = trial.clu(origIdx);
    for kClu = 1 : nClus
        kRes  = res(clu == cluIdx(kClu));
        sc(kClu, :) = histc(kRes, binEdges);
    end
end
