function out = PVBootStrap(gt, varargin)
%   out = PVBootStrap(gt, varargin)
%  [IF_COMPUTE, roi, arena, nResample, nSmplBins, nStd, IF_REPORTFIG]
% pv - population vector , nCells x nSpatialsBins x nTimeBins
    [IF_COMPUTE, IF_REPORTFIG, roi, arena, nResample, nSmplBins, nStd] = ...
        DefaultArgs(varargin, { 0, 0,'CA3', 'bigSquare', 20, 1e3, 5});
    
    if ~IF_COMPUTE
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat']);
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), 'PopVecTimeCourse.mat'], 'dotProd');
        if out.pVal > 0.025, fprintf('\n not sngf \n'); return; end
        figHndl = figure;
        axHdl = subplot(2, 1, 1);
        plot(axHdl, dotProd, 'k');
        axis tight;
        hold on;
        l1Hdl = line(xlim, [out.std_5, out.std_5], 'color', 'r');
        l2Hdl = line(xlim, [out.std_3, out.std_3], 'color', 'g');
        %        legend([l1Hdl, l2Hdl], {'5 std', '3 std'});
        xlabel('# theta cycles', 'FontSize', 14);
        ylabel(' Fisher z - dot product', 'FontSize', 14);
        axHdl = subplot(2, 1, 2);
        logCounts = log10(out.counts);
        logCounts(isinf(logCounts)) = nan;
        bar(out.binEdges, logCounts, 'k');
        hold on;
        axis tight;
        l1Hdl = line([out.std_5, out.std_5], ylim, 'color', 'r');
        l2Hdl = line([out.std_3, out.std_3], ylim, 'color', 'g');
        legend([l1Hdl, l2Hdl], {'5 std', '3 std'});
        xlabel('Fisher z - dot product', 'FontSize', 14);
        ylabel('log_{10}(counts)', 'FontSize', 14, 'interpreter', 'tex');
        if IF_REPORTFIG
            reportfig(figHndl, [mfilename, GenFiletag(roi, arena), gt.datasetType], 0, ['Fb: ', gt.filebase, '<br>', 'trial : ' gt.trialName, '<br>'  'p = ' num2str(out.pVal)]);
        end
        close(figHndl);
        
    else
        %    if isempty(gt.clu), gt = gt.LoadCR;
        if isempty(gt.pfObject), gt = gt.LoadPF; end

        load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
        %    load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), 'PopVecTimeCourse.mat']);
        [popVec, ~, dotProd] = PopVecTimeCourse(gt, [],[],roi, arena);
        [nCells, nXBins, nYBins,  nTimeBins] = size(popVec);
        nSpatialsBins = nXBins * nYBins;
        % popVec = reshape(popVec, nCells, nSpatialsBins, []);
        dpFunc = @(a, b) (a' * b) ./ (norm(a) .* vnorm(b)); % normalized dot product
        smoothedRatMaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus));
        smoothedRatMaps = smoothedRatMaps(:);
        %        oldStr = [];
        matlabpool open local 2
        parfor mResample = 1 : nResample
            pv = full(popVec);
            % scramble cell id & time
            pvs = sparse(pv(randperm(nCells), :, randperm(nTimeBins)));
            dps(:, mResample) = atan(dpFunc(smoothedRatMaps, pvs));
        end
        matlabpool close force local;

        %
        fprintf('\n DONE ..... \n');
        dps(isnan(dps)) = 0;
        binEdges = linspace(0, max(dps(:)), nSmplBins);
        counts = histc(dps(:), binEdges);
        
        sngfVals = nStd * std(dps(:));
        pVal = sum(dps(:) >= sngfVals) / length(dps(:));
        out.binEdges = binEdges; 
        out.counts = counts;
        out.pVal = pVal;
        out.nResample = nResample;
        out.std_5 = sngfVals;
        out.std_3 = 3 * std(dps(:));
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat'], 'out');
    end
end