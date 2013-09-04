function thetaBoundaries = ThetaBoundaries(trial, ThPh,  varargin)
    % thetaBoundaries = ThetaBoundaries(trial, ThPh,  varargin)
    % returns the boundaries of complete theta cycles, boundries identified as the phase with least spiking
    % -------
    % Inputs:
    %    ThPh        - Theta phase time series 
    %    commonClus
    %    thBinSiz
    %    tolerence 
    %    nCatCycles  - number of theta cycles to concatinate
    %    inPeriods   - theta boundaries returned lie with in the inPeriods

    [commonClus, thBinSiz, tolerence, nCatCycles, inPeriods] = DefaultArgs(varargin, {[], 10, 1e-1, 0, []});
    
    res = trial.res(ismember(trial.clu, commonClus)); % load res only for the units in roi
    clu = trial.clu(ismember(trial.clu, commonClus));
    res = round(res .* trial.lfpSampleRate ./ trial.sampleRate) + 1; % convert res to lfp sample rate
    %    [res, resIdx] = SelectPeriods(res, trial.trialPeriods, 'd');
    %  clu = clu(resIdx);
    [res, resIdx] = SelectPeriods(res, trial.trialPeriods , 'd', 1, 1);
    clu = trial.clu(resIdx);
    [trialThPh, thIdx] = SelectPeriods(ThPh, trial.trialPeriods, 'c');
    nBins = 360 / thBinSiz;
    binEdges = linspace(-pi, pi, nBins);
    [count, binIdx] = histc(trialThPh(res), binEdges);
    xx = [binEdges, binEdges(2 : end) + 2*pi];
    [~, minCIdx] = min(count(2:end-1)); 
    minPh = 0.5 * (binEdges(minCIdx) + binEdges(minCIdx + 1));
    thetaBoundaries = find(IsEqual(trialThPh, minPh, tolerence, 0)); 
    if nCatCycles
        thetaBoundaries = thetaBoundaries(1 : nCatCycles : end);
    end
    thetaBoundaries = [thetaBoundaries(1 : end - 1), thetaBoundaries(2 : end)]; 
    if ~isempty(inPeriods)
        thBndInPeriods = IntersectRanges(thetaBoundaries, inPeriods);
        % select only complete theta cycles within inPeriods
        thetaBoundaries = thetaBoundaries(ismember(thetaBoundaries, thBndInPeriods, 'rows'), :);
    end
end