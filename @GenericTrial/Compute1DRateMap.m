function rateMap1D  = Compute1DRateMap(gt, IF_COMPUTE, varargin)
    % rateMap1D  = Compute1DRateMap(gt, IF_COMPUTE, varargin)
    % computes 1D place fields in a linear maze
    % -------
    % Inputs:
    %    IF_COMPUTE  - overwrite
    %    cluId       - cluster ids to consider
    %    nBins       - number of spatial bins
    %    smoothSigma - std of gaussian kernal used to compute smooth PF
    %    posRange    - range of linear positions to consider
    % ------
    % Outputs:
    %    rateMap1D -  cell 1-by-nClus

    if ~IF_COMPUTE & FileExists([gt.paths.analysis,  gt.filebase, '.', gt.trialName, '.', mfilename, '.mat'])
        load([gt.paths.analysis,  gt.filebase, '.', gt.trialName, '.', mfilename, '.mat']);
        return;
    else
        fprintf('\n computing 1D ratemaps... ');
    end
    [res, clu, pos, psp]  = gt.LoadStateRes('RUN', 1, gt.trackingSampleRate);
    [~, linPos] = cart2pol(pos(:, 1), pos(:, 2));
    defPosRange = [min(linPos), max(linPos)];
    defCluId = 1 : NClusters(gt);
    [cluId, nBins, smoothSigma, posRange] = ...
        DefaultArgs(varargin,{defCluId, 50, 3e-2, defPosRange});

    if isempty(gt.clu), gt = gt.LoadCR; end
    binSiz = diff(posRange) ./ nBins;
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