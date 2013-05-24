function y = PVBootStrap(gt, varargin)
%  y = PVBootStrap(pv, dp, varargin
% pv - population vector , nCells x nSpatialsBins x nTimeBins
    [roi, arena, nResample, nSmplBins, nStd, type, IF_REPORTFIG] = ...
        DefaultArgs(varargin, {'CA3', 'bigSquare', 2e1, 1e3, 5, 'display', 0});
    
    if strcmp(type, 'display')
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat']);
        
        if IF_REPORTFIG
            reportfig();
        end
        return;
    end
    %    if isempty(gt.clu), gt = gt.LoadCR;
    if isempty(gt.pfObject), gt = gt.LoadPF; end

    load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
    %    load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), 'PopVecTimeCourse.mat']);
    [popVec, ~, dotProd] = PopVecTimeCourse(gt);
    [nCells, nXBins, nYBins,  nTimeBins] = size(popVec);
    nSpatialsBins = nXBins * nYBins;
    popVec = reshape(popVec, nCells, nSpatialsBins, []);
    dpFunc = @(a, b) (a' * b) ./ (norm(a) .* vnorm(b)); % normalized dot product
    smoothedRatMaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus));
    smoothedRatMaps = smoothedRatMaps(:);
    %        oldStr = [];
    matlabpool open local 2
   parfor mResample = 1 : nResample
       %  str = sprintf(['#', num2str(mResample) ' of ' num2str(nResample)]);
        %   fprintf([repmat('\b', 1, length(oldStr)), str]);
        %o%ldStr = str;
        % scramble in time 
        % pvs  = popVec(:, :, randperm(nTimeBins));
        %scramble cell id
        pvs = popVec(randperm(nCells), :, :);
        dps(:, mResample) = atan(dpFunc(smoothedRatMaps, reshape(pvs, [], nTimeBins )));
        %        dpt(isnan(dpt)) = 0;
        % dps(:, mResample) = dpt;
    end
    matlabpool close;
    %
    dps(isnan(dps)) = 0;
    binEdges = linspace(0, max(dps(:)), nSmplBins);
    counts = histc(dps(:), binEdges);
    
    sngfVals = nStd * std(dps(:));
    pVal = sum(dps(:) >= sngfVals) / length(dps(:));
    
    out.counts = counts;
    out.pVal = pVal;
    out.std_5 = sngfVals;
    out.std_3 = 3 * std(dps(:));
    keyboard;
    save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat'], 'out');
end