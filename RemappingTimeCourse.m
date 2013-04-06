function RemappingTimeCourse(gt, varargin)
    % RemappingTimeCourse(gt, varargin)


    [binSize, roi, tolerence] = DefaultArgs(varargin, { 10,  {'CA1','CA3'}, 1.5e-2});
    %binSize in degrees
    % .1 ~= 5 degrees
    
    load ([gt.paths.data, gt.filebase '.thpar.mat'] );
    if isempty(gt.res)
        gt =gt.LoadCR;
    end
    clus = cell2mat(gt.GetRegionClu(roi)); % shd also add 
    res = gt.res(CompareVectors(gt.clu, clus)); % load res only for the units in roi
    res = round(res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
    [res, resIdx] = SelectPeriods(res, gt.trialPeriods, 'd');
    clu = gt.clu(resIdx);
    [res, resIdx] = SelectPeriods(res, round(gt.goodPosPeriods .* gt.lfpSampleRate ...
        ./ gt.trackingSampleRate) + 1 + res(1),'d');
    clu = clu(resIdx);
    stateThPh = ThPh(res);
    nBins = 360 / binSize;
    bins = linspace(-pi, pi, nBins);
    [count, binIdx] = hist(stateThPh, bins);
    xx = [bins, bins(2:end)+2*pi, bins(end) + 2*pi + 2 * pi / nBins];
    [~, minCIdx] = min(count);
    minPh = xx(minCIdx);
    thetaBoundaries = find(IsEqual(stateThPh, minPh, tolerence)); % in lfp sample rate
    popVec = zeros(length(clus), length(thetaBoundaries) - 1);
    for mClu = 1 : length(clus)
        refVector(mClu) =length(res(clu == clus(mClu))); % divide by time spent
    end
    
    for kPopVec = 1 : length(thetaBoundaries) - 1
        % population vector for each theta cycle, nDims -by- nClycles 
        % nDims = nClus
    
        for kClu = 1 : length(clus)
            curRes = res(clu == clus(kClu));
            popVec(kClu, kPopVec) = sum(curRes >= thetaBoundaries(kPopVec) & curRes <= thetaBoundaries(kPopVec + 1));
        end
    end
    filename = 
    save([gt.paths.analysis, gt.filebase, gt.trialName, 
    %% figures
    bar(xx, [count, count],'FaceColor', 'k');
    hold on;
    line([minPh, minPh], ylim, 'Color', 'g');
    line([minPh, minPh]+2*pi, ylim, 'Color', 'g');

    


end