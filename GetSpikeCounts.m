function [sc, winEdges] = GetSpikeCounts(trial, varargin)
    %  returns spike counts in wins
    % trial - generic trial class object
    % winSize in secs
    % overlap - overlap factor of windows
    % function sc = GetSpikeCounts(trial, winSize, cluIdx)
    
    [winSiz, cluIdx, state, overlap] =  ...
        DefaultArgs(varargin, {200e-3, [], 'RUN', 0});
    %   xySegLength = size(trial.position, 1);
%  winSiz = round(winSiz * trial.sampleRate);
%  nWins = ceil(xySegLength ./ (trial.trackingSampleRate * winSize));
%  winEdges = linspace(0, xySegLength / trial.trackingSampleRate, nWins);
    nClus = length(cluIdx);
    switch trial.datasetType
      case 'MTA'
       %  for i = 1:length(trial.states)
%             if strcmp(state, 'RUN')
%                 state = 'theta';
%             end
%             % if strcmp(trial.states{i}.label, state),
%                 % statePeriods = trial.states{i}.statePeriods;
%             break;
%         end
%         end
       
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
    startEdges = statePeriods(1) : winSiz * (1 - overlap) : statePeriods(2) - winSiz;
    endEdges = statePeriods(1) + winSiz : winSiz * (1 - overlap) : statePeriods(2);
    winEdges = [startEdges', endEdges'] / trial.trackingSampleRate;
    nWins = size(winEdges, 1);
    sc = zeros(nClus, nWins);
    [res, origIdx] = SelectPeriods(trial.res, statePeriods, 'd', 1, 1);
    res = res ./ trial.sampleRate; % res in seconds
    clu = trial.clu(origIdx);
    for mWin = 1 : nWins
        for kClu = 1 : nClus
            kRes  = res(clu == cluIdx(kClu));
            %            sc(kClu, mWin) = histc(kRes, winEdges(mWin, :));
            sc(kClu, mWin) = sum(kRes >= winEdges(mWin, 1)  & kRes <= winEdges(mWin, 2));
        end
    end
end
