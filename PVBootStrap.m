function y = PVBootStrap(gt, varargin)
%  y = PVBootStrap(pv, dp, varargin
% pv - population vector , nCells x nSpatialsBins x nTimeBins
    [roi, arena, nResample] = DefaultArgs(varargin, {'CA3', 'bigSquare', 1e2});



    load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
    %    load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), 'PopVecTimeCourse.mat']);
    [popVec, ~, dotProd] = PopVecTimeCourse(gt);
    [nCells, nXBins, nYBins,  nTimeBins] = size(popVec);
    nSpatialsBins = nXBins * nYBins;
    popVec = reshape(popVec, nCells, nSpatialsBins, []);
    dpFunc = @(a, b) (a' * b) ./ (norm(a) .* vnorm(b)); % normalized dot product
    smoothedRatMaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus));
    smoothedRatMaps = smoothedRatMaps(:);
        oldStr = [];
    %    matlabpool open local 8
    for mResample = 1 : nResample
        str = sprintf(['#', num2str(mResample) ' of ' num2str(nResample)]);
        fprintf([repmat('\b', 1, length(oldStr)), str]);
        oldStr = str;
        % scramble in time 
        pvs  = popVec(:, :, randperm(nTimeBins));
        %scramble cell id
        pvs = popVec(randperm(nCells), :, :);
        dpt = atan(dpFunc(smoothedRatMaps, reshape(pvs, [], nTimeBins )));
        dpt(isnan(dpt)) = 0;
        dps(:, mResample) = dpt;
    end
    
    keyboard;    
end