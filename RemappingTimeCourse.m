

    load ([gt.paths.data, gt.filebase '.thpar.mat'] );
    if isempty(gt.res)
        gt =gt.LoadCR;
    end
    res = round(gt.res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
    res = SelectPeriods(res, gt.trialPeriods, 'd');
    res = SelectPeriods(res, round(gt.goodPosPeriods .* gt.lfpSampleRate ...
        ./ gt.trackingSampleRate) + 1 + res(1),'d');
    stateThPh = ThPh(res);
    binSize = 20; % theta phase in degrees
    nBins = 360 / binSize;
    bins = linspace(-pi, pi, nBins);
    [count, binIdx] = hist(stateThPh, bins);
    
    xx = [bins, bins(2:end)+2*pi, bins(end) + 2*pi + 2 * pi / nBins];
    bar(xx, [count, count],'FaceColor', 'k');
    hold on;
    [~, minCIdx] = min(count);
    minPh = xx(minCIdx);
    line([minPh, minPh], ylim, 'Color', 'g');
    line([minPh, minPh]+2*pi, ylim, 'Color', 'g');
    
    
       
    