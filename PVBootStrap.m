function y = PVBootStrap(gt, pv, dp, varargin)
%  y = PVBootStrap(pv, dp, varargin
% pv - population vector , nCells x nSpatialsBins x nTimeBins
[roi, arena, nResample] = DefaultArgs(varargin, {'CA3', 'bigSquare', 1e1});


   [nCells, nSpatialsBins, nTimeBins] = size(pv);
   load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
   dpFunc = @(a, b) (a' * b) ./ (norm(a) * vnorm(b)); % normalized dot product
   smoothedRatMaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus));
   smoothedRatMaps = smoothedRatMaps(:);
   oldStr = [];
   for mResample = 1 : nResample
       str = sprintf(['#', num2str(mResample) ' of ' num2str(nResample)]);
       fprintf([repmat('\b', 1, length(oldStr)), str]);
       oldStr = str;
       % scramble in time 
       pvs  = pv(:, :, randperm(nTimeBins));
       %scramble cell id
       pvs = pv(randperm(nCells), :, :);
       dps(:, mResample) = arctan(dpFunc(smoothedRatMaps, reshape(pvs, [], nTimeBins)));
   end
   keyboard;    
end