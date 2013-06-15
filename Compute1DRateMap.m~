function rateMap1D  = Compute1DRateMap(gt, varargin)
% rateMap1D  = Compute1DRateMap(gt, res, linPos, varargin)
% res in tracking sample rate

    [res, clu, pos, psp]  = gt.LoadStateRes('RUN', 1, gt.trackingSampleRate);
    [~, linPos] = cart2pol(pos(:, 1), pos(:, 2));
    defPosRange = [min(linPos), max(linPos)];
    if isempty(gt.clu), gt = gt.LoadCR; end
    defCluId = unique(gt.clu);
    [IF_COMPUTE, cluId, nBins, smoothSigma, posRange] = ...
        DefaultArgs(varargin,{0, defCluId, 50, 3e-2, defPosRange});
    if ~IF_COMPUTE
        load([gt.paths.analysis,  gt.filebase, '.', gt.trialName, '.', mfilename, '.mat']);
        return;
    end
    binSiz = diff(defPosRange) ./ nBins;
    edges = linspace(posRange(1), posRange(2), nBins);
    [posCount, binnedPos] = histc(linPos, edges);
    occupancy = posCount ./ gt.trackingSampleRate;
    rateMap1D = cell(length(cluId), 1);
    for kClu = 1 : length(cluId)
        kRes = res(clu == cluId(kClu));
        fprintf('\n cell # %d of %d cells', kClu, length(cluId));
        try
            indx = kRes(find(kRes > 0 & kRes <= length(binnedPos)));
            spkCount = Accumulate([ones(length(indx), 1), binnedPos(indx)], 1, [1, nBins]);
            x = (-nBins : nBins) ./ nBins / 2;
            gw = exp(-x .^ 2 / smoothSigma ^ 2 / 2);
            gw = gw ./ sum(gw);
            smOccupancy = conv(occupancy, gw, 'same');
            smSpkCnt = conv(spkCount, gw, 'same');
            rateMap1D{kClu} = smSpkCnt ./ smOccupancy';
        catch err
            
        end
    end
    save([gt.paths.analysis,  gt.filebase, '.', gt.trialName, '.', mfilename, '.mat'], 'rateMap1D');
end